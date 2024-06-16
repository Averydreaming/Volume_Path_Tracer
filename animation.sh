#!/bin/bash

# cd build
# cmake --build .
# python animation.py

# 设置场景文件夹路径
# scene_folder="./scenes/rayleigh_test/sunrise_to_midday"
scene_folder="$1"

# 确保输出目录存在
output_dir="${scene_folder}/output"
mkdir -p "$output_dir"

# 遍历XML文件
cd "$scene_folder"
for xml_file in *.xml; do
    # 跳过template.xml文件
    if [[ "$xml_file" != "template.xml" ]]; then
        # 构建输出文件名
        echo $xml_file
        x_value=$(echo "$xml_file" | sed -E 's/.*id_([0-9]+)_.*/\1/')
        formatted_x_value=$(printf "%04d" "$x_value")
        output_file="${output_dir}/${formatted_x_value}.exr"
        
        # 构建并执行命令
        cmd="./build/lajolla ${scene_folder}/${xml_file} -o ${output_file}"
        echo "Executing: $cmd"
        $cmd
        
        # 打印进度信息
        progress=$((progress + 1))
        echo "Rendered: $xml_file"
        echo "Finished $progress renders out of $(ls -l "$scene_folder" | grep -v template.xml | wc -l)"
        echo " "
    fi
done

ffmpeg -f image2 -framerate 3 -i "${output_dir}/%04d.exr" -c:v libx264 -crf 18 -pix_fmt yuv420p output.mp4
