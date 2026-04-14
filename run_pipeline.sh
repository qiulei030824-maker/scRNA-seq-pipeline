#!/bin/bash

# ============================================================================
# scRNA-seq Pipeline Runner
# 完整的单细胞分析流程运行脚本
# ============================================================================

set -e  # 遇到错误时退出

# ============================================================================
# 颜色定义
# ============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ============================================================================
# 函数定义
# ============================================================================

print_header() {
    echo -e "${BLUE}"
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo -e "${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

check_command() {
    if command -v $1 &> /dev/null; then
        print_success "$1 已安装"
        return 0
    else
        print_error "$1 未安装"
        return 1
    fi
}

check_file() {
    if [ -f "$1" ]; then
        print_success "文件存在: $1"
        return 0
    else
        print_error "文件不存在: $1"
        return 1
    fi
}

check_directory() {
    if [ -d "$1" ]; then
        print_success "目录存在: $1"
        return 0
    else
        print_warning "目录不存在: $1"
        mkdir -p "$1"
        print_success "已创建目录: $1"
        return 0
    fi
}

# ============================================================================
# 参数解析
# ============================================================================

usage() {
    echo "使用方法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  -i, --input FILE      输入文件 (Seurat .rds 或 10X目录)"
    echo "  -o, --output DIR      输出目录 (默认: ./results)"
    echo "  -s, --species SPECIES 物种 (默认: human)"
    echo "                        可选: human, mouse, arabidopsis, grape"
    echo "  -m, --module MODULE   运行特定模块 (1-9)"
    echo "  -c, --config FILE     配置文件 (默认: config/pipeline_config.yaml)"
    echo "  -h, --help            显示帮助信息"
    echo ""
    echo "示例:"
    echo "  $0 -i data/input.rds -o results/ -s human"
    echo "  $0 -i data/10x_data/ -o results/ -s mouse -m 3"
    echo "  $0 --help"
    exit 1
}

# 默认参数
INPUT_FILE=""
OUTPUT_DIR="./results"
SPECIES="human"
MODULE="all"
CONFIG_FILE="config/pipeline_config.yaml"

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -s|--species)
            SPECIES="$2"
            shift 2
            ;;
        -m|--module)
            MODULE="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_error "未知参数: $1"
            usage
            ;;
    esac
done

# ============================================================================
# 初始化
# ============================================================================

print_header "scRNA-seq Pipeline Runner"
echo "版本: 1.0"
echo "日期: $(date)"
echo ""

# 检查必要参数
if [ -z "$INPUT_FILE" ]; then
    print_error "必须指定输入文件"
    usage
fi

# 检查输入文件
if [ ! -e "$INPUT_FILE" ]; then
    print_error "输入文件不存在: $INPUT_FILE"
    exit 1
fi

# 创建输出目录
check_directory "$OUTPUT_DIR"
check_directory "$OUTPUT_DIR/logs"
check_directory "$OUTPUT_DIR/qc_plots"
check_directory "$OUTPUT_DIR/clustering"
check_directory "$OUTPUT_DIR/annotation"

# 检查配置文件
if [ ! -f "$CONFIG_FILE" ]; then
    print_warning "配置文件不存在: $CONFIG_FILE"
    print_warning "使用默认配置"
else
    print_success "使用配置文件: $CONFIG_FILE"
fi

# ============================================================================
# 检查系统依赖
# ============================================================================

print_header "检查系统依赖"

# 检查R
check_command "R" || {
    print_error "请安装R: https://cran.r-project.org/"
    exit 1
}

# 检查Rscript
check_command "Rscript" || {
    print_error "Rscript不可用"
    exit 1
}

# 检查Python (可选)
if check_command "python3"; then
    PYTHON_AVAILABLE=true
else
    print_warning "Python3未安装，部分功能可能受限"
    PYTHON_AVAILABLE=false
fi

# ============================================================================
# 检查R包依赖
# ============================================================================

print_header "检查R包依赖"

check_r_package() {
    local package=$1
    if Rscript -e "if(!require('$package', quietly=TRUE)) stop('Package $package not found')" &> /dev/null; then
        print_success "R包: $package"
        return 0
    else
        print_error "R包未安装: $package"
        return 1
    fi
}

# 核心R包
CORE_R_PACKAGES=("Seurat" "ggplot2" "dplyr" "patchwork" "Matrix")

for pkg in "${CORE_R_PACKAGES[@]}"; do
    check_r_package "$pkg" || {
        print_warning "尝试安装 $pkg..."
        Rscript -e "install.packages('$pkg', repos='https://cloud.r-project.org')" && print_success "安装成功: $pkg" || print_error "安装失败: $pkg"
    }
done

# ============================================================================
# 运行模块
# ============================================================================

print_header "开始分析流程"
echo "输入文件: $INPUT_FILE"
echo "输出目录: $OUTPUT_DIR"
echo "物种: $SPECIES"
echo "模块: $MODULE"
echo ""

LOG_FILE="$OUTPUT_DIR/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
echo "日志文件: $LOG_FILE"
echo ""

# 记录开始时间
START_TIME=$(date +%s)

# 模块运行函数
run_module() {
    local module_num=$1
    local module_script="scripts/$(printf "%02d" $module_num)_*.R"
    local module_name=""
    
    case $module_num in
        1) module_name="数据加载与质量控制" ;;
        2) module_name="标准化与特征选择" ;;
        3) module_name="降维与聚类" ;;
        4) module_name="标记基因识别" ;;
        5) module_name="AUCell注释" ;;
        6) module_name="PCMaster注释" ;;
        7) module_name="scPlantAnnotate注释" ;;
        8) module_name="集成评分" ;;
        9) module_name="细胞类型命名" ;;
        *) module_name="未知模块" ;;
    esac
    
    print_header "模块 $module_num: $module_name"
    
    # 查找脚本文件
    local script_file=$(find scripts -name "$(printf "%02d" $module_num)_*.R" | head -1)
    
    if [ -z "$script_file" ] || [ ! -f "$script_file" ]; then
        print_error "脚本文件不存在: $script_file"
        return 1
    fi
    
    print_success "运行脚本: $script_file"
    
    # 运行R脚本
    local module_log="$OUTPUT_DIR/logs/module_${module_num}_$(date +%Y%m%d_%H%M%S).log"
    
    if [ $module_num -eq 1 ]; then
        # 模块1需要输入文件和输出目录
        Rscript "$script_file" "$INPUT_FILE" "$OUTPUT_DIR" --species "$SPECIES" 2>&1 | tee "$module_log"
    else
        # 其他模块使用前一个模块的输出
        local prev_output="$OUTPUT_DIR/seurat_$(case $module_num in
            2) echo "filtered" ;;
            3) echo "normalized" ;;
            4) echo "clustered" ;;
            *) echo "processed" ;;
        esac).rds"
        
        if [ ! -f "$prev_output" ]; then
            print_error "前一个模块的输出文件不存在: $prev_output"
            return 1
        fi
        
        Rscript "$script_file" "$prev_output" "$OUTPUT_DIR" 2>&1 | tee "$module_log"
    fi
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ]; then
        print_success "模块 $module_num 完成"
        return 0
    else
        print_error "模块 $module_num 失败"
        return 1
    fi
}

# 运行指定模块或所有模块
if [ "$MODULE" = "all" ]; then
    print_header "运行完整流程 (9个模块)"
    
    for i in {1..9}; do
        if ! run_module $i; then
            print_error "流程在模块 $i 失败"
            exit 1
        fi
        echo ""
    done
    
    print_success "所有模块完成!"
else
    # 运行单个模块
    if [[ "$MODULE" =~ ^[1-9]$ ]]; then
        if ! run_module $MODULE; then
            print_error "模块 $MODULE 运行失败"
            exit 1
        fi
    else
        print_error "无效的模块号: $MODULE (必须是1-9)"
        exit 1
    fi
fi

# ============================================================================
# 生成报告
# ============================================================================

print_header "生成分析报告"

generate_report() {
    local report_script="scripts/generate_report.R"
    
    if [ -f "$report_script" ]; then
        print_success "生成HTML报告..."
        Rscript "$report_script" "$OUTPUT_DIR" 2>&1 | tee "$OUTPUT_DIR/logs/report_generation.log"
        
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            print_success "报告生成完成"
            return 0
        else
            print_warning "报告生成失败"
            return 1
        fi
    else
        print_warning "报告生成脚本不存在: $report_script"
        return 1
    fi
}

generate_report || print_warning "跳过报告生成"

# ============================================================================
# 完成
# ============================================================================

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

print_header "分析完成"
echo "📊 分析总结:"
echo "  输入文件: $INPUT_FILE"
echo "  输出目录: $OUTPUT_DIR"
echo "  物种: $SPECIES"
echo "  运行模块: $MODULE"
echo "  总耗时: $((DURATION / 60))分$((DURATION % 60))秒"
echo ""

echo "📁 输出文件:"
find "$OUTPUT_DIR" -type f -name "*.rds" -o -name "*.csv" -o -name "*.txt" -o -name "*.png" | head -10 | while read file; do
    echo "  - $(basename "$file")"
done

if [ $(find "$OUTPUT_DIR" -type f | wc -l) -gt 10 ]; then
    echo "  ... 和其他 $(($(find "$OUTPUT_DIR" -type f | wc -l) - 10)) 个文件"
fi

echo ""
echo "📈 下一步:"
echo "  1. 检查输出文件: ls -la $OUTPUT_DIR/"
echo "  2. 查看质控图: open $OUTPUT_DIR/qc_plots/*.png"
echo "  3. 检查聚类结果: open $OUTPUT_DIR/clustering/"
echo "  4. 查看注释结果: open $OUTPUT_DIR/annotation/"
echo ""

print_success "scRNA-seq分析流程完成! 🎉"

# 保存运行信息
RUN_INFO="$OUTPUT_DIR/run_info.txt"
cat > "$RUN_INFO" << EOF
scRNA-seq Pipeline Run Information
==================================
运行时间: $(date)
输入文件: $INPUT_FILE
输出目录: $OUTPUT_DIR
物种: $SPECIES
运行模块: $MODULE
总耗时: $((DURATION / 60))分$((DURATION % 60))秒
日志文件: $LOG_FILE
配置文件: $CONFIG_FILE
EOF

print_success "运行信息已保存到: $RUN_INFO"

exit 0