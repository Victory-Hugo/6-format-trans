def replace_zeros_with_n(input_file_path, output_file_path):
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    # 用于存储修改后的内容
    modified_lines = []

    for line in lines:
        line = line.strip()

        if line.startswith(">"):  # 这是标题行，不做修改
            modified_lines.append(line)
        else:
            # 替换序列中的 '0' 为 'N'
            modified_sequence = line.replace('0', 'N')
            modified_lines.append(modified_sequence)

    # 将修改后的内容写入新的文件
    with open(output_file_path, 'w') as output_file:
        for modified_line in modified_lines:
            output_file.write(modified_line + "\n")

    print(f"文件已处理并保存到：{output_file_path}")

# 使用示例
input_file_path = 'C:/Users/LuzHu/Desktop/823_SEA.fasta'  # 请替换为实际的输入文件路径
output_file_path = 'C:/Users/LuzHu/Desktop/823_SEA.N.fasta'  # 请替换为你希望保存的输出文件路径
replace_zeros_with_n(input_file_path, output_file_path)
