#!/bin/bash

# usage: 
# bash run.sh -i <input_file.xml> [-o file_name]
# exp. bash run.sh -i ./scenes/volpath_test/volpath_test1.xml [-o test1] => get a file "test1.exr"
# exp. bash run.sh -i ./scenes/rayleigh_test/sunrise_to_sunset/id_251_sun_120.5_deg.xml

# 设置默认的输入和输出文件名
default_input="input.xml"
default_output="output"

# 初始化输入和输出变量
input_file=$default_input
output_file=$default_output

# 处理命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            if [[ -n "$2" && $2 != -* ]]; then  # 检查是否提供了输入文件且不是另一个选项
                input_file="$2"
                shift # 跳过输入文件参数
            else
                echo "No input file specified for -i/--input"
                exit 1
            fi
            ;;
        -o|--output)
            if [[ -n "$2" && $2 != -* ]]; then  # 检查是否提供了输出文件且不是另一个选项
                output_file="$2"
                shift # 跳过输出文件参数
            else
                echo "No output file specified for -o/--output"
                exit 1
            fi
            ;;
        *) # 未知选项
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift # 跳过当前选项
done

# 执行你的程序，使用输入和输出文件
cmake --build .
./build/lajolla "$input_file" -o "$output_file".exr


# ~/Daily/3-Spring/CG/hdrview/build/HDRView