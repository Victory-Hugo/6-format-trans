import pandas as pd
import numpy as np
from pathlib import Path

def mirror_csv(path_in: str | Path,
               path_out: str | Path,
               name_col: str = 'Unnamed: 0') -> None:
    """
    读取 Fst CSV，将左下三角镜像到右上三角，保留第一列（人群名称），并写回新文件。

    参数
    ----
    path_in  : 输入 CSV 路径
    path_out : 输出 CSV 路径
    name_col : 第一列列名；默认 'Unnamed: 0'
    """
    # ---------- 读取 & 清理 ----------
    df = pd.read_csv(path_in, skipinitialspace=True)          # 去掉列名前的空格
    df[name_col] = df[name_col].str.strip()                   # 行名去空格
    df.columns = df.columns.str.strip()                       # 再保险：列名多余空格
    numeric = df.iloc[:, 1:]                                  # 纯数值区

    # ---------- 检查方阵 ----------
    n_rows, n_cols = numeric.shape
    if n_rows != n_cols:
        raise ValueError(f"除第一列外应为方阵，但得到 {n_rows} × {n_cols}。")

    # ---------- 镜像下三角 ----------
    arr = numeric.to_numpy()                                  # 与 df 同步内存
    triu = np.triu_indices(n_rows, k=1)                       # 严格上三角索引
    arr[triu] = arr.T[triu]                                   # 一行代码完成镜像

    # ---------- 写回 ----------
    df.iloc[:, 1:] = arr
    df.to_csv(path_out, index=False)
    print(f"√ 镜像完成，已保存到 {path_out}")

# --- 快速调用示例 ---
if __name__ == "__main__":
    mirror_csv(
        "/mnt/c/Users/Administrator/Desktop/1.CSV",
        "/mnt/c/Users/Administrator/Desktop/revised_matrix.csv",
    )
