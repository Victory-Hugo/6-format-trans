import pandas as pd
import numpy as np
from pathlib import Path

def full_to_lower_csv(path_in: str | Path,
                      path_out: str | Path,
                      name_col: str = 'Unnamed: 0',
                      keep_diag: bool = False,
                      decimals: int = 7) -> None:
    """
    读取一个全 Fst 矩阵 CSV / TSV，将其裁剪成左下半矩阵（保留第一列），并写回文件。

    参数
    ----
    path_in   : 输入文件路径
    path_out  : 输出文件路径
    name_col  : 第一列列名（默认 'Unnamed: 0'）
    keep_diag : 是否保留对角线(=0)；False 时置为空
    decimals  : 小数保留位数
    """
    # ---------- 读取 & 清理 ----------
    df = pd.read_csv(path_in, sep=r'\s+|,|\t', engine='python', skipinitialspace=True)
    df[name_col] = df[name_col].str.strip()
    df.columns = df.columns.str.strip()
    numeric = df.iloc[:, 1:]

    # ---------- 方阵检查 ----------
    n_rows, n_cols = numeric.shape
    if n_rows != n_cols:
        raise ValueError(f"除第一列外应为方阵，但得到 {n_rows} × {n_cols}。")

    # ---------- 裁剪上三角 ----------
    arr = numeric.to_numpy()
    mask = np.tril(np.ones_like(arr, dtype=bool), k=0 if keep_diag else -1)
    arr[~mask] = np.nan                       # 上三角 & 对角(可选) → NaN
    numeric.iloc[:, :] = arr                  # 写回 DataFrame

    # ---------- 输出 ----------
    float_fmt = lambda x: f"{x:.{decimals}f}" if pd.notna(x) else ''
    df.to_csv(path_out,
              sep='\t',
              index=False,
              header=False,            # 按你的原始需求：无表头
              float_format=float_fmt)

    print(f"√ 半矩阵已保存至: {path_out}")


# -------- 快速调用示例 --------
if __name__ == "__main__":
    full_to_lower_csv(
        "C:/Users/victo/Desktop/新建 Text Document.txt",
        "C:/Users/victo/Desktop/lower_triangle_matrix.txt",
        name_col='Unnamed: 0',   # 如列名并非 'Unnamed: 0' 请自行改成实际列名
        keep_diag=False,         # True 则保留对角 0；False 置空
        decimals=7               # 小数位
    )
