# Orbital_Design_Library

## 概要
このコードは円制限３体問題の軌道設計を行うためのJuliaコードである。

## 動作確認
### 基本情報
- OS : Windows 11 Pro
- CPU : AMD Ryzen Threadripper PRO 5975WX 32-Cores        3.60 GHz
- RAM : 256 GB
- GPU : NVIDIA RTX A6000

### 実行環境に関する詳細
- juliaのバージョン：julia 1.11.3+0
- VScodeのバージョン：

### 動作内容に対する負荷情報
- シミュレーション規模

- 計算時間

- CPU負荷

- GPU負荷

- 消費メモリ

- 

## 計算手法はフォルダ内の以下を参照
Dynamical_Systems_Theory_in_CR3BP.pptx

## 環境構築　(Windowsの場合)
1. Juliaのインストール
Microsoft StoreからJuliaをインストールする。または、[Julia公式ウェブサイト](https://julialang.org/downloads/)からJuliaをダウンロードする。
2. ライブラリのインストール
インストールしたJuliaアプリを実行し、REPLを開く。
3. "]"を入力し、パッケージモードに入り、以下を実行する。
- CUDAのインストール (GPU処理を行うためのライブラリ）
  
```add CUDA```
- CSVのインストール（CSVを扱うライブラリ）
  
```add CSV```
- DataFramesのインストール（データフレーム形式でデータを操作するライブラリ）
  
```add DataFrames```
- DifferentialEquationsのインストール（微分方程式を解くためのライブラリ）
  
```add DifferentialEquations```
- StaticArraysのインストール（固定サイズの配列を扱うためのライブラリ）
  
```add StaticArrays```
- LinearAlgebraのインストール（線形代数計算のライブラリ）
  
```add LinearAlgebra```
- GLMakieのインストール（高性能なデータ可視化ライブラリ）
  
```add GLMakie```
