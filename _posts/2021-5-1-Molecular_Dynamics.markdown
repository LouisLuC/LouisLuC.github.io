---
layout: post
title: "使用GROMACS进行分子动力学模拟案例"
subtitle: 'Example of Moleculer Dynamics'
author: "louis"
mathjax: true
header-style: text
tags:
  - Bioinformatics
  - Moleculer Dynamics
  - CADD
---

## 材料

### GROMACS程序

这里的GROMACS的版本为2021的3月8日的release(2021.1 released March 8th, 2021)

### mdp参数文件

GROMACS的MD中针对不同的步骤如为体系添加离子平衡电荷等需要多套不同的参数，所有参数文件来自于[Bioinformatics-Review](https://github.com/Bioinformatics-Review)的github仓库[MD-simulation-files](https://github.com/Bioinformatics-Review/MD-simulation-files)或者其改变的文件。

使用git将其仓库clone到本地。

```shell
$ git clone https://github.com/Bioinformatics-Review/MD-simulation-files.git
```

本次使用了其中如下的几个文件

- em.mbp
- ions.mbp
- heat.mbp

### 蛋白质PDB文件

本次PSN2的模型使用了自己之前用modeller建立的模型。

将PSN2的模型文件载入pymol中，查看其分子构型

![PSN2_model](../img/in-post/MD_try/PSN2_model.png)

<center>PSN2的模型使用pymol展示。使用多模板建模得到，可以看到其后有个尾巴</center>

这个蛋白质除了中约290-350残基的区域上是一个loop外，其在1-75残基上还有条尾巴，非常长。由于本次模拟仅仅是一次学习，为了加快模拟的速度，防止尾巴和镜像的蛋白相互作用，这次就将其切掉。或者也可以使用更好的模型 (如大名鼎鼎的 [Zhang Lab I-TASSER](https://zhanglab.dcmb.med.umich.edu/I-TASSER/) 建出来的模型)

使用pymol选择1-74个氨基酸，右键delet删除即可

![PSN2_clip](../img/in-post/MD_try/PSN2_clip.png)

切掉后，将其 Export 出来，保存为PSN2_processed.pdb备用。

## 方法

### 产生体系结构和拓扑文件

使用下面的命令从`.pdb`结构文件中产生`.gro`后缀的体系结构和`.top`的拓扑文件。这些文件描述了体系的原子以及它们的位置，还有体系的各个键和电荷等等属性。

```shell
$ gmx pdb2gmx -f PSN2_processed.pdb -o psn2.gro -p psn2.top -ff amber99sb -water tip3p
```
参数的意义如下: 

- -ff: 所选择的立场，立场将用于
- -water: 所选用的水的模型

产生了如下的文件

- psn2.gro
- psn2.top
- posre.itp

如果GROMACS运行成功，那么你从屏幕回显中会得到一个名言警句

> GROMACS reminds you: "We'll celebrate a woman for anything, as long as it's not her talent." (Colleen McCullough)

还有一个信息需要注意

> Total charge -12.000 e

查看输出的屏幕显示，发现电荷为-12 (在`.top`文件的[ atom ]的最后也能看到体系的电荷qtot -12)

### 周期性边界和模拟盒子

由于边界效应会影响独立的系统的边界上粒子的运行性质，MD之前我们得使用周期性边界来消除这种效应。使用如下命令为体系定义边界盒子，这里定义三个盒子，比较它们之中的水分子的量的大小

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

### 加入溶剂

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

### 平衡电荷

在生物体内单独的蛋白理论上应该是一个不带电荷的体系，考虑这一点应该把体系的离子平衡掉

使用下面的命令

```shell
$ gmx grompp -f MD-simulation-files/ions.mdp -c psn2-t-w.gro -p psn2-t.top -o ions.tpr -maxwarn 100
```

- -f: 参数文件
- -c: .gro或者其他结构文件
- -p: .top拓扑文件
- -maxwarn: 设置允许的warnings的文件。默认wen jian允许waranings太低的话将不会输出文件
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
$ gmx genion -s ions.tpr -o psn2-t-w-i.gro -p psn2-t.top -neutral -pname NA -nname CL
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

打开 ***em.mbp*** 文件，建其中的 emtol 更改为500.0，即将停止优化的最大的力的阈值设为 500 *kJ/mol/nm*

再次使用 ***grompp*** 程序，使用参数文件 ***em.mdp*** 产生原子级别的体系描述文件

```bash
$ gmx grompp -f MD-simulation-files/em.mdp -c psn2-t-w-i.gro -p psn2-t.top -o em.tpr
```

根据这套体系描述文件进行能量最小化

```sh
$ gmx mdrun -v -deffnm em
```

- `-v`: 显示每一步的过程以及该过程中体系的一些参数，这些参数包括 
  - `Epot`:  potential energy. should be negative, and (for a simple protein in water) on the order of 10<sup>5</sup>-10<sup>6</sup>
  - `Fmax`: maximum force, Fmax, the target for which was set in em.mdp by the *emtol* argument.
  - `atom`: 原子的数目

跑完后会回显总的能量优化结果，如下

```
Steepest Descents converged to Fmax < 500 in 2153 steps
Potential Energy  = -8.1771481e+05
Maximum force     =  4.8815460e+02 on atom 4500
Norm of force     =  1.0691059e+01
```

将能量最小化过程中的每一步的能量画图展示出来

```shell
$ gmx energy -f em.edr -o potential.xvg
```

执行命令后，程序会要求选择绘制如键、键角以及能量等选项，这里选择10 Potenial 来查看写入能量的数据。

得到的potential.xvg可以使用Xmgrace软件画出来，也可以用编辑器打开，把数据提取出来使用R或者python+matplotlib画出来

本体系选择使用R和ggplot2来画

```R
library(ggplot2)
po = read.table("potential.txt", header = T) # this file is extracted from potential.xvg
p <- ggplot(data=po, mapping=aes(x=Time.ps., y=Potential.kJ.mol.)) + geom_line() + labs(x="potential(kJ/mol)", y="time(ps)")
```



<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/potential.png">    <br>    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系能量最小化时的势能随时间变化的曲线</div> </center>

从上面的曲线图可以看到在能量最小化的过长中体系的能量的确降低收敛了。



### 平衡体系

平衡体系是一个比较难以理解的一部，这一部分涉及的热力学知识超过我的知识水平，有点勉强

首先介3种常用系综: 

- NVT: 正则系综，体系可以和热源交换能量，粒子数固定
- NPT: 巨正则系综，可以交换能量和粒子
- NVE: 能量和离子固定的系统

在我们进行能量最小化后，蛋白质的构象应该变得更为合理，但是如果这时直接开始模拟，体系可能回崩溃。这是因为溶剂分子其分布以及蛋白的取向等等还需要在我们要模拟的温度或者压力下建立合适的取向。

这时候我们可以使用NVT系综来进行的参数文件来构建描述文件，并执行mdrun，由于NVT系综其离子数压力不变，体系是一个升温的过程，并使用一开始得到拓扑文件的 ***posre.itp*** 文件来限制蛋白质中的非氢的原子(也就是所谓的重原子)施加一个除非力很大，都不会移动的限制。

#### NVT平衡

NVT系综中可以和外界热源进行能量的交换，粒子数量固定，我们可以据此来为体系升温。

可以修改MD-simulation-files中的 ***nvt.mdp*** 文件中的参数来控制升温过程中的参数，比较重要的有: 

- `gen_temp` 升温温度(大概?)，默认即为300K，28摄氏度左右，我们的体系就升温到这个温度，不需要更改。
- `nsteps` 执行，基本上决定了体系平衡的时间。计算公式$ 2*x = \frac{1}{1000}time(ps)$ 。默认为50000(100ps)，一般小体系所需时间较少，我们就设置为15000 (即30ps)

注意: 从GROMACS2018开始，参数需要加上 -r 来指定限制性文件，可以和 -c 指定的文件相同(这个参数我不太明白，限制文件不是  ***posre.itp*** 吗?)

> Fatal error:
> Cannot find position restraint file restraint.gro (option -r).
> From GROMACS-2018, you need to specify the position restraint coordinate files
> explicitly to avoid mistakes, although you can still use the same file as you
> specify for the -c option.

```shell
$ gmx grompp -f MD-simulation-files/nvt.mdp -c em.gro -r em.gro -p psn2-t.top -o nvt.tpr
$ gmx mdrun -deffnm nvt
```

输出这样的note

> NOTE 1 [file MD-simulation-files/nvt.mdp]:
> Removing center of mass motion in the presence of position restraints might cause artifacts. When you are using position restraints to equilibrate a macro-molecule, the artifacts are usually negligible.

使用energy模块画出压力和温度变化的数据

```shell
$ gmx energy -f nvt.edr -o temp.xvg
# choose 16 for energy
```





<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/temp.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系升温时的温度随时间变化的曲线</div> 
</center>



体系的温度从一开始就在300K，一直在300K左右震荡，这说明体系的

#### NPT平衡

除了稳定体系的温度，我们还需要稳定体系的压力。使用MD-simulation-files中的npt.mdp文件来进行平衡

这次不更改步数，步数为原始文件中的 50000 ，即100ps

```shell
$ gmx grompp -f MD-simulation-files/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p psn2-t.top -o npt.tpr
$ gmx mdrun -deffnm npt
```

 注意, 我们现在要使用`-t`选项以包括NVT平衡过程中的产生的检查点文件. 这个文件包含了继续模拟所需要的所有状态变量. 为使用NVT过程中得到的速度我们必须包含这个文件. 坐标文件(`-c`)是NVT模拟的最终输出文件.

> NOTE 1 [file psn2-t.top, line 56770]:
>   You are combining position restraints with Parrinello-Rahman pressure
>   coupling, which can lead to instabilities. If you really want to combine
>   position restraints with pressure coupling, we suggest to use Berendsen
>   pressure coupling instead.
>
> NOTE 2 [file MD-simulation-files/npt.mdp]:
>   Removing center of mass motion in the presence of position restraints
>   might cause artifacts. When you are using position restraints to
>   equilibrate a macro-molecule, the artifacts are usually negligible.

查看密度数据和压力的变化

```shell
$ gmx energy -f npt.edr -o density.xvg 
# choose 24 for density
$ gmx energy -f npt.edr -o pressure.xvg 
# choose 18 for pressure
```



<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/pressure.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系压力随时间变化的曲线</div> 
</center>

压力的落差起伏比较大，平均为 $-17.31163 bar$ 。这不能说明有问题，还需要看密度的数据。

<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/density.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系密度随时间变化的曲线</div> 
</center>

密度的数据显示密度平均为 $1015.24 kg/m^3$ ，也在该值附近波动比较稳定，说明我们的体系平衡效果较好，可以同于模拟。

### 运行MD

终于，准备了这么久，我们终于可以进行模拟了🎉🎉

md.mdp文件将其按照上面的公式修改成自己体系需要的时间即可。考虑到性能和时间问题，针对这个体系，将`nsteps`修改为500000，即为1ns;  `dt` 步长设置为0.002，即2fs; 将 `constraints` 更改为 all-bonds 来防止一些原子乱动。

```shell
$ gmx grompp -f MD-simulation-files/md.mdp -c npt.gro -t npt.cpt -p psn2-t.top -o md_1ns.tpr
$ gmx mdrun -v -deffnm md_1ns
```

用`-v`来在屏幕回显上显示执行的步数和预计完成的时间，1ns一般会运行3个小时以上

运行打开vmd或者pymol，选择导入新原子，导入输出的 ***md_1ns.gro*** 文件，之后右键选择导入数据到原子，导入轨迹文件 ***md_1ns.xtc*** 或者 ***md_1ns.trr***，在representaion选项中在展示一栏选择将展示为new cartoon，播放即可看到蛋白质在模拟期间的运动情况。

使用vmd的插件movie maker，就可以将这段轨迹录制下来，比如录制为gif格式

<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/PSN2_MD.gif"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">MD过程中的运动情况.gif</div> 
</center>


从MD的轨迹上来看，除了头尾，其中间的loop的结构变化得最大。

## 分析

MD完成后，我们可以计算MD过程体系的的各个指标以及分析，来探究我们所感兴趣的问题

比较流行的指标有体系原子距离的均方根偏差(RMSD)、均方根涨落(RMSF)、针对氢键的分析、蛋白质的回旋半径的变化以及模拟最终构象与初始构象的差别等等。我们将上面列出来的项目一一分析

### RMSD

使用gmx的rms模块可以直接计算体系模拟过程中的均方根偏差(相对于模拟开始的结构)。使用如下命令指定md前的 ***.tpr*** 文件以及输出的轨迹文件，绘制出C-alpha的RMSD

```shell
$ gmx rms -s md_1ns.tpr -f md_1ns.xtc -o rmsd.xvg -tu ns
# choose 3 for C-alpha for least squares fit
# choose 3 for C-alpha for RMSD calculation
```

其中，`-tu` 时将输出的时间单位定为ns。

使用R画出的图如下

<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/rmsd.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系相对于md开始前的构象的RMSD在MD中随着时间变化的曲线</div> 
</center>
从上图中可以看到，在MD的过程中，RMSD在0.5ns左右即基本稳定在0.4～0.5之间，说明蛋白质的结构在该条件下基本稳定。

### RMSF

输入如下命令，计算体系的C-alpha的均方根涨落。

```shell
$ gmx rmsf -s md_1ns.tpr -f md_1ns.xtc -o rmsf.xvg
# choose 3 for C-alpha
```

<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/rmsf.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系在MD中的均方根涨落</div> 
</center>
除了开头和最后的一些原子外，中间约3500到4700左右的原子的均方根涨落要比周围区域高，这段区域正好对应于 290到350左右的长卷曲loop区域，其涨落更高反映了这个区域的结构松散易动的特性。

### 氢键

氢键由于其数量大，键能不算高而在一定条件下可以断裂的特性，在一些蛋白质的功能的行使中扮演着重要的角色(如配体的结合等等)

现在GROMACS中默认的认定氢键的电子供体和受体的条件如下
$$
r < r_{HB} = 0.35nm \\ \alpha \leqslant \alpha_{HB} = 30
$$

我们使用如下的参数分析整个蛋白质中氢键数量的变化情况

```shell
$ gmx hbond -s md_1ns.tpr -f md_1ns.xtc -num hbond.xvg
# choose 1 and 1 for protein to protein
```

用R将图画出来如下

<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/hbond.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">蛋白质自身的氢键在MD过程中的变化情况</div> 
</center>

可以看到蛋白质在1ns的过程中其氢键数目变化不大，这符合我们MD轨迹观察到的蛋白质除了中间的loop区域之外其他区域构象变化不大的情况。

### 回旋半径

蛋白质的回旋半径可以用于衡量蛋白质的折叠的紧密度。蛋白质的折叠构象足够稳定的话，其回旋半径的变化应该不会大。

使用如下代码得到关于回旋半径的数据，并使用R绘图

```shell
$ gmx gyrate -s md_1ns.tpr -f md_1ns.xtc -o gyrate.xvg
# choose 1 and 1 to calculate hydrogen bonds between protein itself 
```

<center>    
    <img style="border-radius: 0.3125em;"     src="../img/in-post/MD_try/gyrate.png"><br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">体系的回旋半径变化</div> 
</center>

从$R_g$值在模拟的过程中没有巨大的变化，基本保持平稳，说明体系在1ns的份围内其折叠结构处于稳定的状态，没有发生剧烈的构象变化

### 可视化观察构象变化

使用pymol载入PSN2_proceesed.pdb以及最后导出的md_1ns.gro，观察MD前后蛋白构象上的变化

```shell
$ pymol
PyMOL>load PSN2_processed.pdb
PyMOL>load md_1ns.gro
PyMOL>align md_1ns, PSN2_processed
```

将两个结构分别设为两个颜色，我们就可以得到下图的结果

<center>
<figure>
<img src="../img/in-post/MD_try/f_l_comp_2.png" />
<img src="../img/in-post/MD_try/f_l_comp_1.png" />
</figure>
<br/><div style="color:orange; border-bottom: 1px solid #d9d9d9;    display: inline-block;    color: #999;    padding: 2px;">MD前后构象叠加比较的结果</div></center>
可以看到，结构上下面一圈螺旋变化得并不显著，而其 297-357 左右区域内的卷曲由于结构松散而其位置和结构有比较大的变化。



## 结语

以上，对PSN2_HUMAN建模出来的结构进行了模拟，并做了简单的分析。



# 参考

1. http://www.mdtutorials.com/gmx/lysozyme/index.html
2. http://jerkwin.github.io/9999/12/01/GROMACS%E7%A8%8B%E5%BA%8F%E6%96%87%E6%A1%A3/
