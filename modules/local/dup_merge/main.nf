process DUP_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(het_dup_sorted), path(hom_dup_sorted)
    tuple val(meta), path(het_dup_merged), path(hom_dup_merged)
    
    output:
    tuple val(meta), path("*results.bed"), emit: dup_merged_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
        
    file_pairs = [("${het_dup_sorted}", "${het_dup_merged}"), ("${hom_dup_sorted}", "${hom_dup_merged}")]
    
    for sorted, merged in file_pairs:
        bed_DUPs = pd.read_csv(sorted, sep='\\t', header=None)
        bed_merged = pd.read_csv(merged, sep='\\t', header=None)
        
        bed_merged.columns = ["#chr", "start", "end", "program_id"]
        bed_merged["sv_type"] = "DUP"
        bed_merged["DUP_ID"] = bed_merged['#chr'].astype(str) + "_" + bed_merged["start"].astype(str) + "_" + bed_merged['end'].astype(str)
        bed_merged["merged_sv_len"] = bed_merged["end"] - bed_merged["start"] - 1
        
        bed_DUPs["DUP_len"] = bed_DUPs[2] - bed_DUPs[1] - 1
        bed_DUPs["Individual_coord"] = bed_DUPs[3].astype(str) + ":" + bed_DUPs[0].astype(str) + ":" + bed_DUPs[1].astype(str) + "-" + bed_DUPs[2].astype(str)
                
        bed_merged = bed_merged.set_index('DUP_ID')
        df_duplications = pd.DataFrame(columns = ['DUP_ID', 'Program', 'ProgramID', 'Individual_coord','Merged_DUP_LEN', 'Individual_DUP_LEN'])
        
        for ind in bed_merged.index:
            program_id = bed_merged.loc[ind, "program_id"].split(",")
            
            for sv_id in program_id:
                new_row = [ind, sv_id.split(".")[0], sv_id, list(bed_DUPs[bed_DUPs[3] == sv_id]["Individual_coord"])[0], bed_merged["merged_sv_len"][ind], list(bed_DUPs[bed_DUPs[3] == sv_id]["DUP_len"])[0]]
                df_duplications.loc[len(df_duplications.index)] = new_row

        df_duplications["Program_DUP_LEN_TOTAL"] = df_duplications.groupby(["DUP_ID","Program"],as_index = False)["Individual_DUP_LEN"].transform('sum')
        df_duplications["Individual_Overlap"] = round(df_duplications["Individual_DUP_LEN"] / df_duplications["Merged_DUP_LEN"] * 100) 
        df_duplications.loc[df_duplications["Individual_Overlap"] < 0.01, "Individual_Overlap"] = "<0.01"
        
        df_duplications["Program_Overlap_SUM"] = round(df_duplications["Program_DUP_LEN_TOTAL"] / df_duplications["Merged_DUP_LEN"] * 100).clip(upper=100)
        df_duplications.loc[df_duplications["Program_Overlap_SUM"] < 0.01, "Program_Overlap_SUM"] = "<0.01"
        df_duplications["Individual_Overlap_ID"] = df_duplications["ProgramID"] + ":" + df_duplications["Individual_Overlap"].astype("str")
        df_duplications["GT"] = bed_DUPs[7]
        df_duplications["VAF"] = bed_DUPs[6].replace(".", None).astype(float)
                
        n_programs = []
        programs = []
        
        for x in bed_merged["program_id"]: 
            ids_ind = (x.split(","))
            ind_prog = [x.split(".")[0] for x in ids_ind]
            n_programs.append((len(set(ind_prog))))
            programs.append(','.join(set(ind_prog)))
        
        bed_merged["n_programs"] = n_programs
        bed_merged["programs"] = programs
        bed_merged["GT"] = df_duplications.groupby("DUP_ID")["GT"].first()
        bed_merged["mean_VAF"] = df_duplications.groupby("DUP_ID")["VAF"].mean()

        df_duplications["VAF"] = df_duplications["VAF"].fillna(".")
        
        multiple_gt = df_duplications.groupby("DUP_ID")["GT"].nunique()[lambda x: x > 1].index
    
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
        
        for program in list(set(df_duplications["Program"])):
            bed_merged[str(program + "_overlap")] = df_duplications[df_duplications["Program"] == program].groupby("DUP_ID")[["Individual_Overlap_ID"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
            bed_merged[str(program + "_total_overlap")] = df_duplications[df_duplications["Program"] == program].groupby("DUP_ID")[["Program_Overlap_SUM"]].first().astype("str")
            bed_merged[str(program + "_coord")] = df_duplications[df_duplications["Program"] == program].groupby("DUP_ID")[["Individual_coord"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
            bed_merged[str(program + "_GT")] = df_duplications[df_duplications["Program"] == program].groupby("DUP_ID")[["GT"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
            bed_merged[str(program + "_VAF")] = df_duplications[df_duplications["Program"] == program].groupby("DUP_ID")[["VAF"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
        
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
