# jmardr_cfradial

気象庁レーダー毎極座標 GPV の GRIB を netCDF (CF/Radial 規約) に変換する

## 動作環境

-   Python3
    -   numpy
    -   netCDF4
-   wgrib2 (`convert_1hempa.py`のみ)

## 実行ファイル

-   `convert_1hempa.py` 単偏波 GRIB 用
-   `convert_2hempa.py` 二重偏波 GRIB 用

## TODO

-   ドップラー速度：パルス繰り返し周波数の格納

## リンク

-   [配信資料に関する仕様 No.13702 レーダー毎極座標データ（レーダーエコー強度 GPV、ドップラー速度 GPV）](https://www.data.jma.go.jp/add/suishin/shiyou/pdf/no13702)
-   [二重偏波気象レーダー極座標データフォーマット (PPI・RHI 観測共通)](https://www.data.jma.go.jp/developer/past_data/app/2hempa_format_20200914.pdf)
-   [CfRadial Data File Format](https://ral.ucar.edu/projects/titan/docs/radial_formats/CfRadialDoc.pdf)
