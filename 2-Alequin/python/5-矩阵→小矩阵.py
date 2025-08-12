#!/usr/bin/env python3
# /doc: 根据给定的提取名单，从方阵矩阵文件中提取对应的子矩阵并保存，无需临时变量

from pathlib import Path
from typing import List
import pandas as pd


def read_matrix(file_path: Path) -> pd.DataFrame:
    """
    读取制表符分隔的矩阵文件，并以首列作为行索引。

    :param file_path: 矩阵文件路径
    :return: DataFrame，行列索引均为样本ID
    """
    return pd.read_csv(file_path, sep='\t', index_col=0)


def read_id_list(file_path: Path) -> List[str]:
    """
    读取提取名单，每行为一个ID。

    :param file_path: 提取名单文件路径
    :return: ID 列表
    """
    return file_path.read_text(encoding='utf-8').splitlines()


def filter_ids(ids: List[str], df: pd.DataFrame) -> List[str]:
    """
    过滤名单，仅保留同时出现在矩阵行和列中的ID。

    :param ids: 原始ID列表
    :param df: 矩阵 DataFrame
    :return: 有效ID列表
    """
    valid = set(df.index).intersection(df.columns)
    return [i for i in ids if i in valid]


def extract_submatrix(df: pd.DataFrame, ids: List[str]) -> pd.DataFrame:
    """
    从方阵中提取指定ID的子矩阵。

    :param df: 原始方阵 DataFrame
    :param ids: 要提取的ID列表
    :return: 子矩阵 DataFrame
    """
    return df.loc[ids, ids]


def save_matrix(df: pd.DataFrame, out_path: Path) -> None:
    """
    将 DataFrame 保存为制表符分隔文件。

    :param df: 要保存的 DataFrame
    :param out_path: 输出文件路径
    """
    df.to_csv(out_path, sep='\t')


def main():
    # —— 请将以下路径替换成您的实际路径 —— #
    matrix_file = Path(r'Matrix.txt') #todo 矩阵文件，文件包括了行名和列名；制表符分割。
    id_list_file = Path(r'list.txt') #todo 提取名单文件，每行一个ID
    output_file = Path(r'小矩阵.txt') #todo 输出文件路径，保存提取的子矩阵

    # 1. 读取原始矩阵和提取名单
    matrix_df = read_matrix(matrix_file)
    id_list = read_id_list(id_list_file)

    # 2. 过滤出有效的ID
    valid_ids = filter_ids(id_list, matrix_df)

    # 3. 提取子矩阵并保存
    sub_df = extract_submatrix(matrix_df, valid_ids)
    save_matrix(sub_df, output_file)

    # 4. 打印结果
    print(sub_df)


if __name__ == '__main__':
    main()
