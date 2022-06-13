# wd-Al

本仓库的数据和脚本用于复现文章[1]里关于平衡态温稠密铝（密度为2.70 g/cm^3，温度为1 eV）的径向分布函数、静态结构因子和动态结构因子（Dynamic Structure Factor, DSF）。

相关的描述可以参考正文和supplemental material

> [1] Q. Zeng et al., Ab initio validation on the connection between atomistic and hydrodynamic description to unravel the ion dynamics of warm dense matter. [Physical Review Research 3, 033116 (2021)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.033116).

---

- `00.set-train`: 温稠密铝的数据集（密度为2.70 g/cm^3，温度为1 eV），为AIMD的轨迹数据，体系包含32原子，包含24781个构型，可自己自行划分测试集、训练集。

- `01.dp-train`：dpkit的输入文件（基于dpkit1.3.2版本的），对于新版本的dpkit 2.x版本，不建议直接使用，可能会报一些错。

- `02.dpmd`：LAMMPS的输入文件，包含初始弛豫`init.in`和产生NVT系综下的原子轨迹`cal.in`，执行顺序见文件`NOTE`。该步骤会得到径向分布函数`rdf.txt`和原子轨迹文件`dump.*`，后者用于后处理计算密度的时空关联函数

- `03.dsf`：python脚本`test.py`通过LAMMPS输出轨迹文件`02.dpmd/dump/dump.*`计算DSF，由于DSF涉及到不同的波矢取值（k），先后运行`00.generate_run.sh`和`01.sub_all.sh`可以在超算平台（以魔方3的lsf作业系统为例）提交不同波矢的DSF计算，实现某种意义的并行。每个case的python计算时间与体系大小有关，体系为`16x16x16`（包含16384原子）时，计算时间约为1个小时
