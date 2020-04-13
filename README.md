# 基本情報
### コンパイル環境
```
gcc version 9.3.0 (Homebrew GCC 9.3.0)

#基本のコンパイルコマンド
gfortran -fbounds-check ${prg}
```

### 減衰振動関連コード
```
kinddef.f90
oscillation.f90
oscillation_compiler.sh	
```

### Lorenz (1963) モデル関連コード
```
kinddef.f90
lorenz63_prm.f90
lorenz63_cal.f90
lorenz63_compiler.sh
```
### 使用した行列計算ライブラリ
LAPACK(Linear Algebra PACKage)
```
逆行列計算で使用した際のコンパイルコマンド
gfortran -fbounds-check lorenz_*.f90 -o ${prg} -I/usr/local/include -llapack -lblas
```

# 作成情報
- 制作開始日　2020年3月21日

- タグの作成日

1. assimilation_tool_ver1.0(減衰振動のデータ同化手法の実装)　2020年4月1日

### 参考サイト
- http://www.itonwp.sci.u-ryukyu.ac.jp/itokosk.html
- http://www.rcs.arch.t.u-tokyo.ac.jp/kusuhara/tips/linux/fortran.html#sec1
- https://qiita.com/AnchorBlues/items/69c1744de818b5e045ab
