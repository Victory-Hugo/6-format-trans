from Bio import SeqIO
import pandas as pd

# 文件路径
fasta_path = '原始.fasta'
df_path = 'mapping.csv'
output_fasta_path = '重命名.fasta'

# 读取Tip-ID对照表
df_tip = pd.read_csv(df_path)
id_map = dict(zip(df_tip['ID'], df_tip['New_ID']))

# 读取fasta，重命名并保存
records = []
for record in SeqIO.parse(fasta_path, 'fasta'):
    new_id = id_map.get(record.id, record.id)  # 如果无对应则保留原ID
    record.id = new_id
    record.description = ''  # 清除描述
    records.append(record)

SeqIO.write(records, output_fasta_path, 'fasta')
print(f'Renamed FASTA saved to {output_fasta_path}')
