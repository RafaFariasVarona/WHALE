process INS_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(het_ins_sorted), path(hom_ins_sorted)
    
    output:
    tuple val(meta), path("*results.bed"), emit: ins_merged_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
        
    for sorted in ["${het_ins_sorted}", "${hom_ins_sorted}"]:
        df = pd.read_csv(sorted, sep="\\t", header=None)
        df.columns = ["#chr", "start", "end", "program_id", "REF", "ALT", "VAF", "GT"]
        df["INS_ID"] = df['#chr'].astype(str) + "-" + df["start"].astype(str) + "-" + df['end'].astype(str)

        df_ins = pd.DataFrame(columns=['INS_ID', 'Program', 'ProgramID', 'Individual_coord', 'Individual_INS_LEN', 'Merged_INS_LEN'])
        df_ins["INS_ID"] = df['#chr'].astype(str) + "-" + df["start"].astype(str) + "-" + df['end'].astype(str)
        df_ins["Program"] = df["program_id"].apply(lambda x: x.split(".")[0])
        df_ins["ProgramID"] = df["program_id"]
        df_ins["Individual_coord"] = df["program_id"] + ":" + df["#chr"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
        df_ins["Individual_INS_LEN"] = [len(x) - 1 for x in df["ALT"]]
        df_ins["Merged_INS_LEN"] = df_ins.groupby(["INS_ID"], as_index=False)["Individual_INS_LEN"].transform('max')
        df_ins["Program_INS_LEN_TOTAL"] = df_ins.groupby(["INS_ID", "Program"], as_index=False)["Individual_INS_LEN"].transform('sum')
        df_ins["Individual_Overlap"] = round(df_ins["Individual_INS_LEN"] / df_ins["Merged_INS_LEN"] * 100)
        df_ins.loc[df_ins["Individual_Overlap"] < 0.01, "Individual_Overlap"] = "<0.01"
        df_ins["Program_Overlap_SUM"] = round(df_ins["Program_INS_LEN_TOTAL"] / df_ins["Merged_INS_LEN"] * 100).clip(upper=100)
        df_ins.loc[df_ins["Program_Overlap_SUM"] < 0.01, "Program_Overlap_SUM"] = "<0.01"
        df_ins["Individual_Overlap_ID"] = df_ins["ProgramID"] + ":" + df_ins["Individual_Overlap"].astype("str")
        df_ins["GT"] = df["GT"]
        df_ins["VAF"] = df["VAF"]

        merged = df.groupby("INS_ID")["program_id"].apply(lambda x: ','.join(x.unique())).reset_index()    
        merged["sv_type"] = "INS"
        merged["merged_sv_len"] = list(df_ins.groupby("INS_ID")["Merged_INS_LEN"].first())
        
        n_programs = []
        programs = []
        
        for x in merged["program_id"]:
            ids_ind = (x.split(","))
            ind_prog = [x.split(".")[0] for x in ids_ind]
            n_programs.append((len(list(set(ind_prog)))))
            programs.append(','.join(set(ind_prog)))

        merged["n_programs"] = n_programs
        merged["programs"] = programs
        merged["GT"] = list(df_ins.groupby("INS_ID")["GT"].first())
        merged["mean_VAF"] = list(df_ins.groupby("INS_ID")["VAF"].mean())

        multiple_gt = df_ins.groupby("INS_ID")["GT"].nunique()[lambda x: x > 1].index
        
        check_het = 'het' in sorted.lower()
        check_hom = 'hom' in sorted.lower()
        
        def update_gt(row):
            if row["INS_ID"] in multiple_gt:
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
        for program in list(set(df_ins["Program"])):
            program_df = df_ins[df_ins["Program"] == program]
            df_prog_temp = program_df.groupby("INS_ID").agg({
                "Individual_Overlap_ID": lambda x: ','.join(x),
                "Program_Overlap_SUM": 'first',
                "Individual_coord": lambda x: ','.join(x),
                "GT": lambda x: ','.join(x.dropna().astype(str)),
                "VAF": lambda x: ','.join(x.dropna().astype(str)),
            })

            df_prog_temp = df_prog_temp.rename(columns={
                "Individual_Overlap_ID": str(program + "_overlap"),
                "Program_Overlap_SUM": str(program + "_total_overlap"),
                "Individual_coord": str(program + "_coord"),
                "GT": str(program + "_GT"),
                "VAF": str(program + "_VAF")
            })

            if df_prog.empty:
                df_prog = df_prog_temp
            else:
                df_prog = df_prog.join(df_prog_temp, how='outer')
        
        merged = merged.merge(df_prog, on='INS_ID', how='left')

        merged[['#chr', 'start', 'end']] = merged['INS_ID'].str.split('-', expand=True)
        
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
