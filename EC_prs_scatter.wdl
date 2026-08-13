version 1.0

workflow plink_prs {

  input {
    Array[String] chromosomes
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    File cll_sumstats_file

    String marker_col_name
    String beta_col_name
    String score_allele_col_name
  }

  scatter (i in range(length(chromosomes))) {

    call RunPRSBatch {
      input:
        chr_name = chromosomes[i],
        pgen = pgen_files[i],
        pvar = pvar_files[i],
        psam = psam_files[i],
        sumstats_file = cll_sumstats_file,
        marker_col_name = marker_col_name,
        beta_col_name = beta_col_name,
        score_allele_col_name = score_allele_col_name
    }
  }

  output {
    Array[File] prs_archives = RunPRSBatch.prs_archive
  }
}

task RunPRSBatch {
  input {
    String chr_name
    File pgen
    File pvar
    File psam
    File sumstats_file
    String marker_col_name
    String beta_col_name
    String score_allele_col_name
  }

  command <<<
    set -euo pipefail

    mkdir -p results/~{chr_name}
    prefix=$(basename "~{sumstats_file}" .tsv)
    outdir="results/~{chr_name}/${prefix}"
    mkdir -p "$outdir"

    header=$(head -n1 ~{sumstats_file})
    marker_col=$(echo "$header" | tr '\t' '\n' | nl -v1 | awk '$2=="~{marker_col_name}"{print $1}')
    allele_col=$(echo "$header" | tr '\t' '\n' | nl -v1 | awk '$2=="~{score_allele_col_name}"{print $1}')
    beta_col=$(echo "$header" | tr '\t' '\n' | nl -v1 | awk '$2=="~{beta_col_name}"{print $1}')

    if [ -z "$marker_col" ] || [ -z "$allele_col" ] || [ -z "$beta_col" ]; then
      echo "ERROR: Failed to detect required columns."
      echo "$header"
      exit 1
    fi

    echo -e "SNP\tA1\tBETA" > "$outdir/${prefix}.weights.tsv"
    awk -v mc="$marker_col" -v ac="$allele_col" -v bc="$beta_col" 'BEGIN{FS=OFS="\t"} NR>1 {print $mc, $ac, $bc}' \
      ~{sumstats_file} >> "$outdir/${prefix}.weights.tsv"

    awk 'NR>1 {print $1}' "$outdir/${prefix}.weights.tsv" > "$outdir/${prefix}.snps_to_extract.txt"

    mkdir -p pfile_input
    ln -s ~{pgen} pfile_input/chr_input.pgen
    ln -s ~{pvar} pfile_input/chr_input.pvar
    ln -s ~{psam} pfile_input/chr_input.psam

    plink2 \
      --pfile pfile_input/chr_input \
      --extract "$outdir/${prefix}.snps_to_extract.txt" \
      --score "$outdir/${prefix}.weights.tsv" 1 2 3 header cols=scoresums \
      --out "$outdir/${prefix}_PRS"

    tar -czf prs_results_~{chr_name}.tar.gz -C results ~{chr_name}
  >>>

  output {
    File prs_archive = "prs_results_~{chr_name}.tar.gz"
  }

  runtime {
    docker: "emosyne/plink2:7.2_Alpha5.9"
    memory: "256 GB"
    cpu: 8
    disks: "local-disk 100 HDD"
  }
}