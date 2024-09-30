process BND_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(het_bnd_sorted), path(hom_bnd_sorted)
    
    output:
    tuple val(meta), path("*results.bed"), emit: bnd_merged_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
        
    for sorted in ["${het_bnd_sorted}", "${hom_bnd_sorted}"]:
        df = pd.read_csv(sorted, sep="\\t", header=None)
        df.columns = ["#chr", "start", "end", "program_id", "REF", "ALT", "VAF", "GT"]
        df["end"] = df["start"] + 1
        df["BND_ID"] = df['#chr'].astype(str) + "-" + df["start"].astype(str) + "-" + df["end"].astype(str)
        df["chr2"] = df['ALT'].str.extract(r'[\\[\\]](.*?):')

        df_bnd = pd.DataFrame(columns=['BND_ID', 'Program', 'ProgramID', 'Individual_coord', 'Merged_BND_LEN', 'Individual_BND_LEN'])
        df_bnd["BND_ID"] = df['#chr'].astype(str) + "-" + df["start"].astype(str) + "-" + df["end"].astype(str)
        df_bnd["chr2"] = df['ALT'].str.extract(r'[\\[\\]](.*?):')
        
        df_bnd["Program"] = df["program_id"].apply(lambda x: x.split(".")[0])
        df_bnd["ProgramID"] = df["program_id"]
        df_bnd["Individual_coord"] = df["program_id"] + ":" + df["#chr"].astype(str) + ":" + (df["start"] + 1).astype(str) + "-" + df['ALT']
        df_bnd["Individual_BND_LEN"] = None
        df_bnd["Merged_BND_LEN"] = None
        df_bnd["Program_BND_LEN_TOTAL"] = None
        df_bnd["Individual_Overlap"] = None
        df_bnd["Program_Overlap_SUM"] = None
        df_bnd["Individual_Overlap_ID"] = None
        df_bnd["GT"] = df["GT"]
        df_bnd["VAF"] = df["VAF"].replace(".", None).astype(float)

        merged = df.groupby(["BND_ID", "chr2"])["program_id"].apply(lambda x: ','.join(x.unique())).reset_index()
        merged["sv_type"] = "BND"
        merged["merged_sv_len"] = list(df_bnd.groupby(["BND_ID", "chr2"])["Merged_BND_LEN"].first())

        n_programs = []
        programs = []

        for x in merged["program_id"]:
            ids_ind = (x.split(","))
            ind_prog = [x.split(".")[0] for x in ids_ind]
            n_programs.append((len(list(set(ind_prog)))))
            programs.append(','.join(set(ind_prog)))

        merged["n_programs"] = n_programs
        merged["programs"] = programs
        merged["GT"] = list(df_bnd.groupby(["BND_ID", "chr2"])["GT"].first())
        merged["mean_VAF"] = list(df_bnd.groupby(["BND_ID", "chr2"])["VAF"].mean())
        
        multiple_gt = df_bnd.groupby(["BND_ID", "chr2"])["GT"].nunique()[lambda x: x > 1].index
        
        check_het = 'het' in sorted.lower()
        check_hom = 'hom' in sorted.lower()
        
        def update_gt(row):
            if (row["BND_ID"], row["chr2"]) in multiple_gt:
                if row["GT"] == "0/0":
                    return "0/1"
                elif row["GT"] == "./.":
                    if check_het:
                        return "0/1"
                    elif check_hom:
                        return "1/1"
            return row["GT"]
        
        merged["GT"] = merged.apply(update_gt, axis=1)
        
        df_prog = pd.DataFrame()
        for program in list(set(df_bnd["Program"])):
            program_df = df_bnd[df_bnd["Program"] == program]
            df_prog_temp = program_df.groupby(["BND_ID", "chr2"]).agg({
                "Individual_coord": lambda x: ','.join(x),
                "GT": lambda x: ','.join(x.dropna().astype(str)),
                "VAF": lambda x: ','.join(x.dropna().astype(str)),
            })

            df_prog_temp = df_prog_temp.rename(columns={
                "Individual_coord": str(program + "_coord"),
                "GT": str(program + "_GT"),
                "VAF": str(program + "_VAF")
            })

            if df_prog.empty:
                df_prog = df_prog_temp
            else:
                df_prog = df_prog.join(df_prog_temp, how='outer')
            
            df_prog[str(program + "_overlap")] = None
            df_prog[str(program + "_total_overlap")] = None
        
        merged = merged.merge(df_prog, on='BND_ID', how='left')

        merged[['#chr', 'start', 'end']] = merged['BND_ID'].str.split('-', expand=True)
    
        cols = ['#chr', 'start', 'end'] + [col for col in merged.columns if col not in ['#chr', 'start', 'end']]
        merged = merged[cols]

        column_order = [
            "#chr", "start", "end", "program_id", "sv_type", "merged_sv_len",
            "n_programs", "programs", "GT", "mean_VAF",
            "svim_overlap", "svim_total_overlap", "svim_coord", "svim_GT", "svim_VAF",
            "cuteSV_overlap", "cuteSV_total_overlap", "cuteSV_coord", "cuteSV_GT", "cuteSV_VAF",
            "Sniffles2_overlap", "Sniffles2_total_overlap", "Sniffles2_coord", "Sniffles2_GT", "Sniffles2_VAF"
        ]

        for col in column_order:
            if col not in merged.columns:
                merged[col] = pd.NA
        
        merged = merged[column_order]

        output_file = os.path.splitext(sorted)[0] + "_merged_results.bed"
        merged.to_csv(output_file, index=False, sep="\\t")
    """
}
