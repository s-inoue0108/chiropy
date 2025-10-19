# CHIROpy

Gaussian TD-DFT 計算のファイル生成とキラル光学特性解析を行うユーティリティプログラムです。

## 依存関係

Python3 と、以下の Python パッケージを必要とします：

- NumPy
- RDKit
- OpenBabel
  - eigen
- NetworkX
- Matplotlib

また、OpenBabel 3.1.\* のバイナリプログラムを別途必要とします。MacOSX であれば、HomeBrew (`brew install openbabel`) で入れることができます。

## インストール

Python 3.10 以上の環境を用意し、wheel バイナリから `pip install` によって依存パッケージごとインストールすることができます。

```bash
pip install chiropy-***.whl
```

```bash
# 確認
pip list | grep chiropy

# CLI のテスト
chiropy --help
```

# Gaussian インプットの生成 (`chiropy ginp`)

1分子の化学構造ファイルから、TD-DFT 計算用の Gaussian インプットを生成します。

```bash
# 化学構造ファイルからインプットを生成する
chiropy ginp -i *.mol

# ヘルプ
chiropy ginp --help
```

化学構造ソースとして利用可能なファイル形式は以下の通りです。

| 拡張子 | 説明 | 備考 |
| :---- | :-- | :--|
| `mol` | MDL mol ファイル (V2000/V3000) | 1分子のもの |
| `sdf` | MDL sdf ファイル (V2000/V3000) | mol ファイルと同等/1分子のもの |
| `mol2` | mol2 ファイル | 1分子のもの |
| `xyz` | xyz ファイル | 1分子のもの |
| `out` | Gaussian アウトプットファイル | 構造最適化後のアウトプットが適する |
| `log` | Gaussian アウトプットファイル | 構造最適化後のアウトプットが適する |

## 励起状態計算

TD-DFT 計算を用いた励起状態計算を行うためのインプットを生成します。

```bash
# 第1励起状態をターゲット・第3励起状態までを探索
# td(root=1, nstate=3) (デフォルト)
chiropy ginp -i *.mol 

# td(root=1, nstate=8)
chiropy ginp -i *.mol --nstate 8

# td(root=2, nstate=4)
chiropy ginp -i *.mol --root 2 --nstate 4
```

## 励起状態の構造最適化

`--opt` オプションを指定すると、`--root` で指定したターゲットの励起状態を構造最適化します。

```bash
# 第1励起状態を最適化
# td(root=1, nstate=3) opt
chiropy ginp -i *.mol --opt

# td(root=1, nstate=8) opt
chiropy ginp -i *.mol --nstate 8 --opt

# td(root=2, nstate=4) opt
chiropy ginp -i *.mol --root 2 --nstate 4 --opt
```

## Gaussian の条件

### DFT の汎関数

`-f` または `--functional` で指定可能です。デフォルトは `b3lyp` です。

```bash
# CAM-B3LYP に変更
nicspy ginp -i *.mol -f "cam-b3lyp"
```

### 基底関数系

`-b` または `--basis` で指定可能です。デフォルトは `6-31g(d)` です。

```bash
# def2-TZVP に変更
nicspy ginp -i *.mol -b "def2tzvp"

# wB97X-D/6-311G(d,p) に変更
nicspy ginp -i *.mol -f "wb97xd" -b "6-311g(d,p)"
```

### 有効内殻ポテンシャル

`--ecp` で指定可能です。指定した ECP を重金属に自動で割り当てます。

```bash
# LanL2DZ を指定
nicspy ginp -i *.mol --ecp "lanl2dz"

# def2-SVP/def2-TZVP を指定
nicspy ginp -i *.mol -b "def2svp" --ecp "def2tzvp"
```

### その他のオプション

- 総電荷の変更: `-c [charge]` (デフォルトは `-c 0`)
- スピン多重度の変更: `-m [multiplicity]` (デフォルトは `-m 1`)
- SCF 収束条件の変更: `--scf tight` など
- 出力ファイル名を指定 `-o [filename]`
- 出力ファイルの接尾辞の変更 `-s [suffix]` (デフォルトは `-s "_tddft"`)
- 出力ファイルの拡張子の変更: `-e com` (デフォルトは `-e gjf`)
- Gaussian 出力のレベルの変更: `--verbose 1` (デフォルトは `--verbose 0`)

# Gaussian アウトプットの解析 (`chiropy gout`)

アウトプットファイルを読み取り、キラル光学特性に関わる物理量を計算します。特に、ETDM および MTDM についてはベクトルを化学構造上にマッピング表示します。

```bash
# アウトプットファイルを解析する (*.out または *.log)
chiropy gout -i *.out

# ヘルプ
chiropy gout --help
```

## GUI による表示

3次元化学構造をレンダリングし、ETDM および MTDM ベクトルを描画して GUI を表示します。

```bash
# GUI で表示
chiropy gout -i *.out
```

### 励起状態の選択

`--state` を指定することで、励起状態を選択することができます。デフォルトは `1` です。

```bash
第3励起状態の情報を表示
chiropy gout -i *.out --state 3
```

### 表示オプション

- `--symbol`: Ball-and-stick モデルの代わりに元素記号を表示 (より軽量にレンダリングできます)
- `--showh`: 水素を表示
- `--notext`: 各種テキストを非表示
- `--dark`: ダークテーマで表示

## その他オプション

- `--image`: GUI を表示せず、画像をエクスポート
- `--summary`: GUI を表示せず、結果の概要を標準出力に表示

