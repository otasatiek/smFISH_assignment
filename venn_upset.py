import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
from matplotlib_venn import venn2, venn3
import argparse
import pathlib

def plot_and_count_intersections(csv_path, output_path, venn_threshold, bin_threshold, out_csv_path=None):
    """
    CSVデータを読み込み、共発現パターンのカウントCSV出力と、
    カテゴリ数に応じたベン図/UpSetプロットの生成を行う。

    Args:
        --csv (str): 入力CSVファイルのパス。
        --output (str): 出力画像ファイルのパス。
        --venn_threshold (int): この数未満のカテゴリ数であればベン図を描画する閾値。
        --bin_threshold (float): データを二値化する際の閾値。この値以上をTrueとする。
        --out_csv (str, optional): カウント結果を保存するCSVファイルのパス。
    """
    try:
        print(f"Loading CSV: {csv_path}")
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"エラー: CSVファイルが見つかりません: {csv_path}")
        return

    # 'nucleus_label' を除くすべての列を処理対象とする
    if 'nucleus_label' not in df.columns:
        print("警告: 'nucleus_label' 列が見つかりません。すべての列を対象とします。")
        target_cols = sorted(df.columns.tolist())
    else:
        target_cols = sorted([col for col in df.columns if col != 'nucleus_label'])
    
    if not target_cols:
        print("エラー: 'nucleus_label' 以外の処理対象列が見つかりませんでした。")
        return

    num_sets = len(target_cols)
    print(f"Found {num_sets} sets for plotting: {', '.join(target_cols)}")
    print(f"Binarization threshold for sets: value >= {bin_threshold}")

    # --- データをbool型（True/False）に変換 ---
    df_bool = pd.DataFrame()
    for col in target_cols:
        df_bool[col] = df[col] >= bin_threshold

    # --- カウントCSVの生成 ---
    if out_csv_path:
        print("Generating co-expression pattern counts...")
        try:
            # groupbyで各パターンの組み合わせの数をカウント
            pattern_counts = df_bool.groupby(target_cols).size().reset_index(name='count')
            
            # 0のカウントは不要な場合が多いため除外
            pattern_counts = pattern_counts[pattern_counts['count'] > 0]

            # カウント数で降順にソート
            pattern_counts = pattern_counts.sort_values('count', ascending=False)
            
            # CSVとして保存
            pattern_counts.to_csv(out_csv_path, index=False)
            print(f"  -> Saved counts to {out_csv_path}")
        except Exception as e:
            print(f"  -> エラー: CSVカウントファイルの保存中に問題が発生しました。詳細: {e}")

    plt.rcParams['font.size'] = 18  # デフォルトのフォントサイズを16に設定
    # --- プロットの分岐ロジック ---
    plt.figure()

    if 1 < num_sets < venn_threshold:
        # **ベン図を生成**
        print(f"Number of sets ({num_sets}) is less than threshold ({venn_threshold}). Generating Venn Diagram.")
        set_names = [col.replace('_spot_count', '').replace('_positive', '') for col in target_cols]
        
        if num_sets == 2:
            set1_data = df_bool[target_cols[0]]
            set2_data = df_bool[target_cols[1]]
            subset_counts = (
                (set1_data & ~set2_data).sum(),
                (~set1_data & set2_data).sum(),
                (set1_data & set2_data).sum()
            )
            venn2(subsets=subset_counts, set_labels=set_names)
            plt.title(f"Intersection of 2 Sets (Threshold >= {bin_threshold})")

        elif num_sets == 3:
            set1_data = df_bool[target_cols[0]]
            set2_data = df_bool[target_cols[1]]
            set3_data = df_bool[target_cols[2]]
            subset_counts = (
                (set1_data & ~set2_data & ~set3_data).sum(),
                (~set1_data & set2_data & ~set3_data).sum(),
                (set1_data & set2_data & ~set3_data).sum(),
                (~set1_data & ~set2_data & set3_data).sum(),
                (set1_data & ~set2_data & set3_data).sum(),
                (~set1_data & set2_data & set3_data).sum(),
                (set1_data & set2_data & set3_data).sum()
            )
            venn3(subsets=subset_counts, set_labels=set_names)
            plt.title(f"Intersection of 3 Sets (Threshold >= {bin_threshold})")
        
    else:
        # **UpSetプロットを生成**
        print(f"Number of sets ({num_sets}) is not less than threshold ({venn_threshold}). Generating UpSet Plot.")
        
        df_upset = df_bool.set_index(pd.Index(target_cols)).transpose()
        plot_data = UpSet(df_upset, subset_size='count', show_counts=True)
        plot_data.plot()
        plt.suptitle(f"Intersection of Sets (Threshold >= {bin_threshold})")

    # プロットをファイルに保存
    try:
        plt.savefig(output_path, bbox_inches='tight')
        print("-" * 40)
        print(f"✅ Plot saved successfully to: {output_path}")
    except Exception as e:
        print(f"エラー: プロットの保存中に問題が発生しました。詳細: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="カテゴリの数に応じて、ベン図またはUpSetプロットと、共発現パターンのカウントCSVを生成します。"
    )
    parser.add_argument('--csv', required=True, help="入力となるデータCSVファイルのパス。'nucleus_label'以外の全数値列が対象です。")
    parser.add_argument('--output', default='plot.png', help='出力画像ファイルのパス (デフォルト: plot.png)')
    parser.add_argument(
        '--out_csv',
        default=None,
        help='(オプション) 共発現パターンのカウントを保存するCSVファイルのパス。指定しない場合は保存されません。'
    )
    parser.add_argument(
        '--venn_threshold',
        type=int,
        default=4,
        help='この数未満のカテゴリ数であればベン図を生成する閾値 (デフォルト: 4)。'
    )
    parser.add_argument(
        '--bin_threshold',
        type=float,
        default=2.0,
        help='集合を作成する際の二値化の閾値。この値以上を陽性(True)とします。(デフォルト: 2.0)'
    )

    args = parser.parse_args()
    
    plot_and_count_intersections(args.csv, args.output, args.venn_threshold, args.bin_threshold, args.out_csv)