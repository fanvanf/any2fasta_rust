#!/usr/bin/env bash
set -euo pipefail

# ===================== 配置项 =====================
PERL_BIN="/home/f/miniforge3/bin/any2fasta"
RUST_BIN="./target/release/any2fasta"
TMP_DIR="./test_tmp"
TEST_CASES=(
    "fasta_plain"
    "fasta_gz"
    "fastq_plain"
    "fastq_bz2"
    "genbank_plain"
    "embl_plain"
    "gfa_plain"
    "clustal_plain"
    "pdb_plain"
)
PARAM_SETS=(
    ""
    "-u"
    "-l"
    "-n -u"
    "-s -u"
    "-g -n"
    "-p -n"
)
# =================================================

# 初始化环境
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"
PASS_COUNT=0
FAIL_COUNT=0

# 用于存储性能数据的关联数组
declare -A PERL_TIME
declare -A PERL_MEM
declare -A RUST_TIME
declare -A RUST_MEM

# 保存所有测试用例名称（按执行顺序）
declare -a TEST_NAMES

info() {
    echo -e "\033[34m[INFO] $1\033[0m"
}
pass() {
    echo -e "\033[32m[PASS] $1\033[0m"
    ((PASS_COUNT++))
}
fail() {
    echo -e "\033[31m[FAIL] $1\033[0m"
    ((FAIL_COUNT++))
}

# 1. 校验程序是否存在
[[ -x "${PERL_BIN}" ]] || { echo "ERROR: Perl脚本不存在: ${PERL_BIN}"; exit 1; }
[[ -x "${RUST_BIN}" ]] || { echo "ERROR: Rust程序不存在，请先编译: cargo build --release"; exit 1; }

# 检查 /usr/bin/time 是否存在
TIME_CMD="/usr/bin/time"
if [[ ! -x "${TIME_CMD}" ]]; then
    echo "WARNING: /usr/bin/time 不存在，跳过性能测试"
    TIME_CMD=""
fi

# 2. 生成各类测试文件
generate_test_files() {
    info "开始生成测试用例文件..."

    # FASTA 核酸
    cat > "${TMP_DIR}/fasta_plain.fasta" <<EOF
>seq1 chromosome 1 geneA
ATGCGGATTCGGNN--ATGC
>seq2 geneB
TTGACCGXAT
EOF
    gzip -c "${TMP_DIR}/fasta_plain.fasta" > "${TMP_DIR}/fasta_gz.fasta.gz"

    # FASTQ
    cat > "${TMP_DIR}/fastq_plain.fq" <<EOF
@read1 illumina:1
AGCTTGCN
+
AAAAA###
@read2
TTGGCC
+
FFFFFF
EOF
    bzip2 -c "${TMP_DIR}/fastq_plain.fq" > "${TMP_DIR}/fastq_bz2.fq.bz2"

    # GenBank
    cat > "${TMP_DIR}/genbank_plain.gbk" <<EOF
LOCUS       TEST001 100 bp DNA linear CON 01-JAN-2025
DEFINITION  test sequence.
VERSION     TEST001.1
ORIGIN
        1 atggcctagg gcttagccgt nnnn--ggcc
//
EOF

    # EMBL/UniProt蛋白
    cat > "${TMP_DIR}/embl_plain.embl" <<EOF
ID   PROT01; SV 1; linear; genomic DNA; STD; UNC; 30 AA.
SQ   SEQUENCE   30 AA;  3500 MW;
     ALA ARG ASN ASP CYS GLY
//
EOF

    # GFA
    cat > "${TMP_DIR}/gfa_plain.gfa" <<EOF
H	VN:Z:1.0
S	11	ACCTT
S	12	TCAAGG
S	13	CTTGATT
L	11	+	12	-	4M
L	12	-	13	+	5M
L	11	+	13	+	3M
P	14	11+,12-,13+	4M,5M
EOF

    # Clustal
    cat > "${TMP_DIR}/clustal_plain.clw" <<EOF
CLUSTAL multiple sequence alignment

seq1   ATGCGG--NN
seq2   TTAGCCXXAT
EOF

    # PDB
    cat > "${TMP_DIR}/pdb_plain.pdb" <<EOF
HEADER TEST PROTEIN 01-JAN-25 1ABC
SEQRES   1 A  6  ALA ARG ASN ASP CYS GLY
SEQRES   1 B  4  LEU ILE VAL MET
EOF

    info "所有测试文件生成完毕"
}

# 3. 执行对比（含性能测量）
run_compare() {
    local test_name="$1"
    local input_file="$2"
    local params="$3"
    local out_perl="${TMP_DIR}/${test_name}_perl.fasta"
    local out_rust="${TMP_DIR}/${test_name}_rust.fasta"
    local log_perl="${out_perl%.fasta}.log"
    local log_rust="${out_rust%.fasta}.log"
    local time_perl="${TMP_DIR}/${test_name}_perl.time"
    local time_rust="${TMP_DIR}/${test_name}_rust.time"
    local md5_perl md5_rust
    local perl_elapsed perl_mem rust_elapsed rust_mem

    info "测试用例: ${test_name} | 参数: [${params}]"

    # ---------- 执行 Perl 版本（计时） ----------
    if [[ -n "${TIME_CMD}" ]]; then
        { time ${PERL_BIN} ${params} "${input_file}" > "${out_perl}" 2> "${log_perl}" ; } 2> "${time_perl}"
        # 解析时间文件（格式: elapsed_seconds max_rss_kb）
        read perl_elapsed perl_mem < "${time_perl}" || true
        PERL_TIME["${test_name}"]="${perl_elapsed:-0}"
        PERL_MEM["${test_name}"]="${perl_mem:-0}"
    else
        # 如果不支持time，直接执行
        ${PERL_BIN} ${params} "${input_file}" > "${out_perl}" 2> "${log_perl}"
        PERL_TIME["${test_name}"]="N/A"
        PERL_MEM["${test_name}"]="N/A"
    fi

    # 检查Perl执行是否成功（如果失败则直接标记失败并跳过Rust）
    if [[ ! -s "${out_perl}" ]] && grep -q "ERROR" "${log_perl}" 2>/dev/null; then
        fail "Perl 执行失败: ${test_name}"
        return 1
    fi

    # ---------- 执行 Rust 版本（计时） ----------
    if [[ -n "${TIME_CMD}" ]]; then
        { time ${RUST_BIN} ${params} "${input_file}" > "${out_rust}" 2> "${log_rust}" ; } 2> "${time_rust}"
        read rust_elapsed rust_mem < "${time_rust}" || true
        RUST_TIME["${test_name}"]="${rust_elapsed:-0}"
        RUST_MEM["${test_name}"]="${rust_mem:-0}"
    else
        ${RUST_BIN} ${params} "${input_file}" > "${out_rust}" 2> "${log_rust}"
        RUST_TIME["${test_name}"]="N/A"
        RUST_MEM["${test_name}"]="N/A"
    fi

    if [[ ! -s "${out_rust}" ]] && grep -q "ERROR" "${log_rust}" 2>/dev/null; then
        fail "Rust 执行失败: ${test_name}"
        return 1
    fi

    # ---------- MD5 比对 ----------
    md5_perl=$(md5sum "${out_perl}" 2>/dev/null | awk '{print $1}')
    md5_rust=$(md5sum "${out_rust}" 2>/dev/null | awk '{print $1}')

    if [[ -z "${md5_perl}" || -z "${md5_rust}" ]]; then
        fail "MD5 计算失败: ${test_name}"
        return 1
    fi

    if [[ "${md5_perl}" == "${md5_rust}" ]]; then
        pass "${test_name} MD5一致"
    else
        fail "${test_name} MD5不一致"
        diff -u "${out_perl}" "${out_rust}" > "${TMP_DIR}/${test_name}_diff.log" 2>/dev/null || true
        info "差异已保存至: ${TMP_DIR}/${test_name}_diff.log"
    fi

    # 记录测试名称用于性能汇总
    TEST_NAMES+=("${test_name}")
    return 0
}

# ========== 主流程 ==========
generate_test_files

# 遍历所有文件+参数组合
for case in "${TEST_CASES[@]}"; do
    case_file=$(find "${TMP_DIR}" -name "${case}.*" -type f ! -name "*.log" ! -name "*.time" ! -name "*.diff" | head -1)
    [[ -z "${case_file}" ]] && continue

    for param in "${PARAM_SETS[@]}"; do
        test_id="${case}_$(echo ${param} | tr ' ' '_')"
        run_compare "${test_id}" "${case_file}" "${param}" || true
    done
done

# ========== 汇总报告 ==========
echo -e "\n==================== 测试汇总 ===================="
echo "总通过: ${PASS_COUNT}"
echo "总失败: ${FAIL_COUNT}"
if [[ ${FAIL_COUNT} -eq 0 ]]; then
    echo -e "\033[32m全部测试用例输出完全一致！\033[0m"
else
    echo -e "\033[31m存在输出不一致，请检查diff日志\033[0m"
fi

# 性能测试结果（如果可用）
if [[ -n "${TIME_CMD}" ]] && [[ ${#TEST_NAMES[@]} -gt 0 ]]; then
    echo -e "\n==================== 性能测试结果 ===================="
    printf "%-45s %15s %15s %15s %15s\n" "Test Case" "Perl Time(s)" "Rust Time(s)" "Perl Mem(KB)" "Rust Mem(KB)"
    for name in "${TEST_NAMES[@]}"; do
        perl_t="${PERL_TIME[$name]:-N/A}"
        rust_t="${RUST_TIME[$name]:-N/A}"
        perl_m="${PERL_MEM[$name]:-N/A}"
        rust_m="${RUST_MEM[$name]:-N/A}"
        printf "%-45s %15s %15s %15s %15s\n" "$name" "$perl_t" "$rust_t" "$perl_m" "$rust_m"
    done

    # 计算并显示加速比（仅当都有数值）
    echo -e "\n加速比 (Rust/Perl 耗时):"
    for name in "${TEST_NAMES[@]}"; do
        pt="${PERL_TIME[$name]}"
        rt="${RUST_TIME[$name]}"
        if [[ "$pt" != "N/A" && "$rt" != "N/A" && "$pt" != "0" && "$rt" != "0" ]]; then
            ratio=$(echo "scale=2; $pt / $rt" | bc -l 2>/dev/null || echo "N/A")
            printf "%-45s %s\n" "$name" "$ratio"
        fi
    done
else
    echo -e "\n性能测试未启用（/usr/bin/time 不存在）"
fi

echo -e "\nDone."

if [[ ${FAIL_COUNT} -eq 0 ]]; then
    exit 0
else
    exit 1
fi