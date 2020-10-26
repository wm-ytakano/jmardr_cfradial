# jmardr_cfradial

気象庁レーダー毎極座標 GPV の GRIB を netCDF (CF/Radial 規約) に変換する

## 動作環境

-   Python3
    -   numpy
    -   netCDF4
-   wgrib2

## TODO

-   偏波パラメータなど他要素への対応
-   volume_number に何を格納すべきか
-   エコー強度：エコーなしが 0dBZ 扱いなのをどう扱うべきか
-   ドップラー速度：パルス繰り返し周波数の格納

## リンク

-   [配信資料に関する仕様 No.13702 レーダー毎極座標データ（レーダーエコー強度 GPV、ドップラー速度 GPV）](https://www.data.jma.go.jp/add/suishin/shiyou/pdf/no13702)
-   [CfRadial Data File Format](https://ral.ucar.edu/projects/titan/docs/radial_formats/CfRadialDoc.pdf)
