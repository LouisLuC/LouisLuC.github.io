## Introduce

这里的GROMACS的版本为2021的3月8日的release(2021.1 released March 8th, 2021)

## preperation

### 准备mdp立场参数

GROMACS的MD中针对不同的步骤如为体系添加离子平衡电荷等需要多套不同的立场参数，所有立场参数文件来自于[Bioinformatics-Review](https://github.com/Bioinformatics-Review)的github仓库[MD-simulation-files](https://github.com/Bioinformatics-Review/MD-simulation-files)，使用了其中的[em.mbp](https://github.com/Bioinformatics-Review/MD-simulation-files/blob/master/em.mdp)以及[ions.mdp](https://github.com/Bioinformatics-Review/MD-simulation-files/blob/master/ions.mdp)。

使用git将其仓库clone到本地。

```shell
$ git clone https://github.com/Bioinformatics-Review/MD-simulation-files.git
```

本次使用了其中如下的几个文件

- em.mbp
- ions.mbp
- heat.mbp

### 准备PDB文件

将PSN2的模型文件载入pymol中，查看其分子构型

![PSN2_model](/home/louislu/repository/LouisLuC.github.io/img/in-post/MD_try/PSN2_model.png)

<center>PSN2的模型使用pymol展示。使用多模板建模得到，可以看到其后有个尾巴</center>

这条尾巴有大约70个氨基酸，非常长。由于本次模拟仅仅是一次学习，为了加快模拟的速度，防止尾巴和镜像的蛋白相互作用，这次就将其切掉。或者也可以使用更好的模型 (如大名鼎鼎的 [Zhang Lab I-TASSER](https://zhanglab.dcmb.med.umich.edu/I-TASSER/) 建出来的模型)

使用pymol选择1-74个氨基酸，右键delet删除即可

![PSN2_clip](/home/louislu/repository/LouisLuC.github.io/img/in-post/MD_try/PSN2_clip.png)

切掉后，将其 Export 出来，保存为PSN2_processed.pdb备用。

## Methods

### use command *gmx pdb2gmx* to get the topology from the pdb file

use command line fellowed to generate topology file of PSN2_model.pdb. 

```shell
$ gmx pdb2gmx -f PSN2_processed.pdb -o psn2.gro -p psn2.top -ff amber99sb -water tip3p
```
the parameters

- -ff: Force field. (option) You can choose is behind
- -water: water

this step will geneterat 3 files: 

- psn2.gro
- psn2.top
- posre.itp

if gromacs command run successfully, you shall see one saying on your screen.

> GROMACS reminds you: "We'll celebrate a woman for anything, as long as it's not her talent." (Colleen McCullough)

Then we retrive three files: ***prosre.itp***, ***psn2.gro*** and ***psn2.top***
the ***prosre.itp*** file is used to definded the position registion regular.

> Total charge -12.000 e

查看输出的屏幕显示，发现电荷为-12 (在.top文件的[ atom ]的最后也能看到体系的电荷qtot -12)

### Create Simulator Box

由于边界效应，MD都使用周期性边界()来消除这种效应。使用如下命令为体系定义边界盒子

```shell
$ gmx editconf -f psn2.gro -o psn2-c.gro -bt cubic -d 0.8 -c
$ gmx editconf -f psn2.gro -o psn2-t.gro -bt triclinic -d 0.8 -c #done
$ gmx editconf -f psn2.gro -o psn2-d.gro -bt dodecahedron -d 0.8 -c
```

- -f: 输入的体系结构文件
- -bt: 加的盒子类型
- -d: 蛋白距离盒子边界的最小距离()
- -c: 设置蛋白质处于中心位置
- -o: 输出文件

为了接下来各个体系之间不相互影响，将top文件分为3份，按盒子类型分为 _psn2-t.top_、_psn2-c.top_、_psn2-d.top_

### Add water solution and minimalize the energy

由于细胞环境为溶液环境，模拟需要添加水环境，使用如下命令为3个体系添加水分子

从这一步开始，每一步将会修改 **_.top_** 文件，可以做好备份，或者程序每一步都会生成一个修改之前的备份 **#XXX.top.n#** 文件，需要记住每一个备份和其步骤的对应关系

```shell
$ gmx solvate -cp psn2-c.gro -cs spc216.gro -p psn2-c.top -o psn2-c-w.gro
$ gmx solvate -cp psn2-t.gro -cs spc216.gro -p psn2-t.top -o psn2-t-w.gro
$ gmx solvate -cp psn2-d.gro -cs spc216.gro -p psn2-d.top -o psn2-d-w.gro
```

注意备份一下.top文件，每次添加分子都会修改.top文件(或者每次命令也会备份一个用"#"框起来的.top文件)

```shell
$ gmx solvate -cp psn2-tri.gro -p psn2-t.top -cs spc216.gro -o psn2-tri-w.gro # done
$ gmx solvate -cp psn2-cubic.gro -p psn2-c.top -cs spc216.gro -o psn2-cubic-w.gro
$ gmx solvate -cp psn2-dodecahedron.gro -p psn2-d.top -cs spc216.gro -o psn2-dodecahedron.gro
```

- -cp: 添加了模拟盒子后的体系文件
- -p: .top拓扑文件
- -cs: 添加的溶剂的结构
- -o: 输出文件

每次添加分子后，屏幕回显会显示添加溶剂分子的量

- d: Number of solvent molecules:  19723
- c: Number of solvent molecules:  28847
- t: Number of solvent molecules:  16253

这些信息也可以在相应的**_.top_**文件中找到

```vim
[ molecules ]
; Compound        #mols
Protein_chain_A     1
SOL             16253
```

为了节省模拟的时间，接下来使用 triclinic 盒子进行模拟

## 平衡电荷

在生物体内单独的蛋白理论上应该是一个不带电荷的体系，考虑这一点应该把体系的离子平衡掉

使用下面的命令

```shell
$ gmx grompp -f MD-simulation-files/ions.mdp -c psn2-t-w.gro -p psn2-t.top -o ions.tpr -maxwarn 100
```

- -f: 参数文件
- -c: .gro或者其他结构文件
- -p: .top拓扑文件
- -maxwarn: 设置允许的warnings的文件。默认允许waranings太低的话将不会输出文件
- -o: .tpr输出，原子级别的体系的描述

这个命令会有WARNING或者NOTE，告知用户该体系的电荷没有屁股很高，以及算法上的注意事项等等

> NOTE 1 [file MD-simulation-files/ions.mdp]:
>   With Verlet lists the optimal nstlist is >= 10, with GPUs >= 20. Note
>   that with the Verlet scheme, nstlist has no effect on the accuracy of
>   your simulation.
>
> Generating 1-4 interactions: fudge = 0.5
>
> NOTE 2 [file psn2-t.top, line 56769]:
>   System has non-zero total charge: -12.000000
>   Total charge should normally be an integer. See
>   http://www.gromacs.org/Documentation/Floating_Point_Arithmetic
>   for discussion on how close it should be to an integer.
>
> 
>
> Number of degrees of freedom in T-Coupling group rest is 115485.00
>
> NOTE 3 [file MD-simulation-files/ions.mdp]:
>   You are using a plain Coulomb cut-off, which might produce artifacts.
>   You might want to consider using PME electrostatics.

从上面的提示可以看到，我们的体系中有12个负电荷，所以要添加12个正电荷来平衡

使用如下的命令为体系添加电荷

```shell
$ gmx genion -s ions.tpr -o psn2-t-w-i.gro -p psn2-t.top -neutral -pname NA
```

- -s: **_.tpr_** 文件
- -neutral: 设置了参数后程序将会自动添加平衡的电荷
- -cc: 添加的离子的浓度
- -pame: 正电荷离子名，不指定的话会添加Na<sup>+</sup> (不是 Not avaliable! )
- -nname: 负电荷离子名

添加完离子后，可以在 **_.top_** 文件中的***[ molecules ]***栏目中看到添加的离子的信息

```
[ molecules ]
; Compound        #mols
Protein_chain_A     1
SOL         16241
NA               12
```

### 能量最小化

加了溶剂以及平衡了体系的电荷之后，在我们进行模拟之前，需要确保体系的几何构型的正确、没有原子之间没有碰撞等问题，放映为不合理的构象以及偏高的能量。经过能量最小化，体系的构象会更加合理。

再次使用 ***grompp***程序产生原子级别的体系描述文件

```bash
$ gmx grompp -f MD-simulation-files/em.mdp -c psn2-t-w-i.gro -p psn2.top -o em.tpr
```



蛋白出盒子

氢键分析选一个区域来算氢键，必须选某个区域所有的原子来算氢键

编辑mdp文件来看信息

官网上找做采样和升温的mdp

水最少的盒子跑

