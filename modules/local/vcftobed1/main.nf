process VCFTOBED1 {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk 'BEGIN { FS="\\t"; OFS="\\t" }
    {
        if (\$0 !~ /^#/)
        {
            info=\$8
            split(info, info_fields, ";")
            end=""
            af=""
            for (i in info_fields)
            {
                split(info_fields[i], kv, "=")
                if (kv[1] == "END")
                {
                    end=kv[2]
                }
                if (kv[1] == "AF")
                {
                    af=kv[2]
                }
            }

            format=\$9
            sample=\$10
            split(format, format_fields, ":")
            split(sample, sample_fields, ":")

            for (i in format_fields)
            {
                if (format_fields[i] == "GT")
                {
                    gt=sample_fields[i]
                }
            }

            print \$1, \$2 - 1, end, \$3, \$4, \$5, af, gt
        }
    }' ${vcf} > ${vcf.getBaseName()}.bed
    """
}
