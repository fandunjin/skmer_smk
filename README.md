# Skmer with snakemake workflow
## 首先，拉取本Repositories
```bash
git@github.com:fandunjin/skmer_smk.git
```
## 接着，创建skmer工作环境
```bash
conda create -n skmer snakemake Jellyfish Mash seqtk pandas=1.5.2 scipy biopython
# skmer correct中pandas的版本需要<2，否则会报错rep91 ValueError，所以需要单独安装较低版本的pandas
conda activate skmer #进入环境
pip show pandas # 确认pandas版本
```
### 由于conda中的skmer版本过低，需要从github中下载：https://github.com/shahab-sarmashghi/Skmer
```bash
git clone https://github.com/shahab-sarmashghi/Skmer.git
cd Skmer
python setup.py install
```
## 最后，将参考基因组数据ref_fna放入raw_data中，原始双端测序数据放入raw_data/raw_data中，后运行snakemake即可。
### 需要注意的是：
#### 1.命名格式要一致：参考基因组为ref_fna；原始双端测序数据为sample.R1.fq.gz和sample.R2.fq.gz（sample可修改，后缀一定要为".R1.fq.gz"）
#### 2.数据放在正确的位置：ref_fna在raw_data中，sample.R1.fq.gz在raw_data/raw_data中
```bash
snakemake --core 48 --latency-wait 120 #回到snakefile所处路径，运行此命令行即可完成分析
```

## 输出结果包含：

## 测试流程
```bash
git@github.com:fandunjin/skmer_smk.git
conda create -n skmer snakemake Jellyfish Mash seqtk pandas=1.5.2 scipy biopython
conda activate skmer
git clone https://github.com/shahab-sarmashghi/Skmer.git
cd Skmer
python setup.py install
cd ../
snakemake --core 4
```