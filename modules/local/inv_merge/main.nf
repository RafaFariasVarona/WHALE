process INV_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(het_inv_sorted), path(hom_inv_sorted)
    tuple val(meta), path(het_inv_merged), path(hom_inv_merged)
    
    output:
    tuple val(meta), path("*results.bed"), emit: inv_merged_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    
    file_pairs = [("${het_inv_sorted}", "${het_inv_merged}"), ("${hom_inv_sorted}", "${hom_inv_merged}")]
    
    for sorted, merged in file_pairs:
        bed_INVs = pd.read_csv(sorted, sep='\\t', header=None)
        bed_merged = pd.read_csv(merged, sep='\\t', header=None)

        bed_merged.columns = ["#chr", "start", "end", "program_id"]
        bed_merged["sv_type"] = "INV"
        bed_merged["INV_ID"] = bed_merged['#chr'].astype(str) + "_" + bed_merged["start"].astype(str) + "_" + bed_merged['end'].astype(str)
        bed_merged["merged_sv_len"] = bed_merged["end"] - bed_merged["start"] - 1
        
        bed_INVs["inv_len"] = bed_INVs[2] - bed_INVs[1] - 1
        bed_INVs["Individual_coord"] = bed_INVs[3].astype(str) + ":" + bed_INVs[0].astype(str) + ":" + bed_INVs[1].astype(str) + "-" + bed_INVs[2].astype(str)
                
        bed_merged = bed_merged.set_index('INV_ID')
        df_inversions = pd.DataFrame(columns = ['INV_ID', 'Program', 'ProgramID', 'Individual_coord','Merged_INV_LEN', 'Individual_INV_LEN'])
        
        for ind in bed_merged.index:
            program_id = bed_merged.loc[ind, "program_id"].split(",")
            
            for sv_id in program_id:
                new_row = [ind, sv_id.split(".")[0], sv_id, list(bed_INVs[bed_INVs[3] == sv_id]["Individual_coord"])[0], bed_merged["merged_sv_len"][ind], list(bed_INVs[bed_INVs[3] == sv_id]["inv_len"])[0]]
                df_inversions.loc[len(df_inversions.index)] = new_row
        
        df_inversions["Program_INV_LEN_TOTAL"] = df_inversions.groupby(["INV_ID","Program"],as_index = False)["Individual_INV_LEN"].transform('sum')
        df_inversions["Individual_Overlap"] = round(df_inversions["Individual_INV_LEN"] / df_inversions["Merged_INV_LEN"] * 100) 
        df_inversions.loc[df_inversions["Individual_Overlap"] < 0.01, "Individual_Overlap"] = "<0.01"
        
        df_inversions["Program_Overlap_SUM"] = round(df_inversions["Program_INV_LEN_TOTAL"] / df_inversions["Merged_INV_LEN"] * 100).clip(upper=100)
        df_inversions.loc[df_inversions["Program_Overlap_SUM"] < 0.01, "Program_Overlap_SUM"] = "<0.01"
        df_inversions["Individual_Overlap_ID"] = df_inversions["ProgramID"] + ":" + df_inversions["Individual_Overlap"].astype("str")
        df_inversions["GT"] = bed_INVs[7]
        df_inversions["VAF"] = bed_INVs[6]
                
        n_programs = []
        programs = []
        
        for x in bed_merged["program_id"]:
            ids_ind = (x.split(","))
            ind_prog = [x.split(".")[0] for x in ids_ind]
            n_programs.append((len(set(ind_prog))))
            programs.append(','.join(set(ind_prog)))
        
        bed_merged["n_programs"] = n_programs
        bed_merged["programs"] = programs
        bed_merged["GT"] = df_inversions.groupby("INV_ID")["GT"].first()
        bed_merged["mean_VAF"] = df_inversions.groupby("INV_ID")["VAF"].mean()
        
        multiple_gt = df_inversions.groupby("INV_ID")["GT"].nunique()[lambda x: x > 1].index
    
        check_het = 'het' in sorted.lower()
        check_hom = 'hom' in sorted.lower()
        
        def update_gt(row):
            if row.name in multiple_gt:
                if row["GT"] == "0/0":
                    return "0/1"
                elif row["GT"] == "./.":
                    if check_het:
                        return "0/1"
                    elif check_hom:
                        return "1/1"
            return row["GT"]
        
        bed_merged["GT"] = bed_merged.apply(update_gt, axis=1)
        
        for program in list(set(df_inversions["Program"])):
            bed_merged[str(program + "_overlap")] = df_inversions[df_inversions["Program"] == program].groupby("INV_ID")[["Individual_Overlap_ID"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
            bed_merged[str(program + "_total_overlap")] = df_inversions[df_inversions["Program"] == program].groupby("INV_ID")[["Program_Overlap_SUM"]].first().astype("str")
            bed_merged[str(program + "_coord")] = df_inversions[df_inversions["Program"] == program].groupby("INV_ID")[["Individual_coord"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
            bed_merged[str(program + "_GT")] = df_inversions[df_inversions["Program"] == program].groupby("INV_ID")[["GT"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
            bed_merged[str(program + "_VAF")] = df_inversions[df_inversions["Program"] == program].groupby("INV_ID")[["VAF"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
        
        column_order = [
            "#chr", "start", "end", "program_id", "sv_type", "merged_sv_len",
            "n_programs", "programs", "GT", "mean_VAF",
            "svim_overlap", "svim_total_overlap", "svim_coord", "svim_GT", "svim_VAF",
            "cuteSV_overlap", "cuteSV_total_overlap", "cuteSV_coord", "cuteSV_GT", "cuteSV_VAF",
            "Sniffles2_overlap", "Sniffles2_total_overlap", "Sniffles2_coord", "Sniffles2_GT", "Sniffles2_VAF"
        ]
        
        for col in column_order:
            if col not in bed_merged.columns:
                bed_merged[col] = pd.NA
        
        bed_merged = bed_merged[column_order]

        output_file = os.path.splitext(sorted)[0] + "_merged_results.bed"
        bed_merged.to_csv(output_file, index=False, sep="\\t")
    """
}
