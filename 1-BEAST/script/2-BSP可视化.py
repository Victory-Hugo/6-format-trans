import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from aquarel import load_theme

theme = load_theme('boxy_light')
theme.apply()

def read_data(file_path: str, smooth: bool):
    """
    读取文件内容，并返回标题、表头和数据。
    
    若 smooth 为 True：保留浮点精度；否则将数据转换为整数。
    """
    content = Path(file_path).read_text(encoding='utf-8').splitlines()
    if len(content) < 3:
        raise ValueError("文件内容不符合预期格式。")
    
    title = content[0].strip()
    header_line = content[1].strip()
    header = header_line.split('\t')
    data_lines = content[2:]
    
    def is_numeric(line: str) -> bool:
        try:
            [float(x) for x in line.strip().split('\t')]
            return True
        except ValueError:
            return False

    # 删除末尾非数字的行
    while data_lines and not is_numeric(data_lines[-1]):
        data_lines.pop()
    
    if smooth:
        # 保持浮点精度
        data = [
            [float(item) for item in line.strip().split('\t')]
            for line in data_lines
        ]
    else:
        # 先转换为float再转换为int
        data = [
            [int(float(item)) for item in line.strip().split('\t')]
            for line in data_lines
        ]
    return title, header, data

def apply_smoothing(df: pd.DataFrame, window_size: int = 3) -> pd.DataFrame:
    """
    应用移动平均平滑，窗口大小默认设为3。
    """
    return df.rolling(
        window=window_size,
        min_periods=1,
        center=True
    ).mean()

def main(bsp_file: str, smooth_str: str, time_lower_limit: int, time_upper_limit: int, OUTPUT_FILE: str):
    """
    参数：
        bsp_file: 输入文件路径
        smooth_str: "YES" 表示采用平滑；"NO" 表示禁用平滑
        time_lower_limit: 时间下限（整数）
        time_upper_limit: 时间上限（整数）
    """
    # 根据传入参数确定是否采用平滑
    smooth_function = True if smooth_str.upper() == "YES" else False

    title, header, data = read_data(bsp_file, smooth_function)
    
    # 创建DataFrame，确保数据为float类型，并以'time'列作为索引
    df = pd.DataFrame(data, columns=header).astype(float)
    df.set_index('time', inplace=True)
    
    if smooth_function:
        df_processed = apply_smoothing(df)
    else:
        df_processed = df.round()
    
    # 根据传入的时间下限和上限进行过滤
    data_filtered = df_processed.loc[time_lower_limit:time_upper_limit]
    
    # 根据模式转换数据类型
    if smooth_function:
        time_values = data_filtered.index.to_numpy()
        upper_values = data_filtered['upper'].to_numpy()
        lower_values = data_filtered['lower'].to_numpy()
    else:
        time_values = data_filtered.index.to_numpy().astype(int)
        upper_values = data_filtered['upper'].to_numpy().astype(int)
        lower_values = data_filtered['lower'].to_numpy().astype(int)
    
    # 计算Y轴范围
    """
    y_min = 10 ** floor(log10(min(lower_values)))
    y_max = 10 ** ceil(log10(max(upper_values)))
    """
    y_min = 10 ** np.floor(np.log10(lower_values.min()))
    y_max = 10 ** np.ceil(np.log10(upper_values.max()))
    
    # 设置字体和生成图形文本为矢量可编辑
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    # 绘图
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(time_values, data_filtered['mean'].round().astype(int),
            label='Mean', color='#C36B5F', linewidth=1, linestyle='--')
    ax.plot(time_values, data_filtered['median'].round().astype(int),
            label='Median', color='#72B6AB', linewidth=1)
    ax.fill_between(time_values, upper_values, lower_values,
                    color='#B9DEDD', alpha=0.3)
    
    ax.grid(False)
    ax.set_title('Bayesian skyline')
    ax.set_xlabel('Time')
    ax.set_ylabel('Ne')
    ax.set_yscale('log')
    ax.set_ylim(y_min, y_max)
    
    # 设置Y轴刻度为10的幂次
    yticks = [10 ** i for i in range(int(np.log10(y_min)), int(np.log10(y_max)) + 1)]
    ax.set_yticks(yticks)
    ax.set_yticklabels([f'{tick:.0f}' for tick in yticks])
    
    ax.legend()
    theme.apply_transforms()
    plt.savefig(OUTPUT_FILE)
    plt.show()

# 示例调用（请根据实际情况传入参数）
if __name__ == '__main__':
    BSP_FILE = "/mnt/c/Users/victo/Desktop/EA.BSP"   # 替换为实际文件路径
    SMOOTH = "YES"                           # "YES" 或 "NO"
    TIME_LOWER_LIMIT = 0                    # 时间下限
    TIME_UPPER_LIMIT = 25000                # 时间上限
    OUTPUT_FILE = "/mnt/c/Users/victo/Desktop/BSP.pdf"

    main(BSP_FILE, SMOOTH, TIME_LOWER_LIMIT, TIME_UPPER_LIMIT,OUTPUT_FILE)