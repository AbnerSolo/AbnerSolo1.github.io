# GAUSSIAN 阅读笔记：从平滑参数到最坏情况—平均情况归约

[返回笔记目录](../../README.md)

整理：AbnerSolo · 2026-10-08  
论文：Daniele Micciancio and Oded Regev, *Worst-case to Average-case Reductions based on Gaussian Measures*。初版发表于 FOCS 2004，期刊版发表于 *SIAM Journal on Computing* 37(1):267–302, 2007。

> **版本与阅读范围。** 本文的定义、引理、定理编号和页码均对应作者公开的 **2005-12-14、33 页版本**：[原文 PDF][paper]。期刊信息见[作者论文主页][homepage]与 [DOI][doi]。本文是学习笔记，不是论文翻译；重点是从基础概念推出平滑性质，再解释它怎样支撑归约。本文中的 GAUSSIAN 不指 GPV 的 *How to Use a Short Basis*（ePrint 2007/432）。
>
> **标记约定。** “严格结论”给出明确假设与结论；“直观解释”帮助建立图像；“证明思路”不替代完整证明；“本文推导”指从已列公式自行推得的补充；“待核实／未展开”列出仍需回到原文处理的范围。核心公式已对照上述作者稿，核对范围见第 10 节。

## 目录

1. [论文要解决什么问题](#sec-1)
2. [符号与基础对象](#sec-2)
3. [对偶格、Fourier 与 Poisson 各负责什么](#sec-3)
4. [平滑参数：定义、意义与精确结论](#sec-4)
5. [怎样用格的几何量控制平滑参数](#sec-5)
6. [为什么还需要离散高斯的矩估计](#sec-6)
7. [高斯如何生成归约需要的样本](#sec-7)
8. [从这些工具到最坏情况—平均情况归约](#sec-8)
9. [论文贡献与容易混淆的说法](#sec-9)
10. [核对记录与待进一步核实的部分](#sec-10)
11. [参考文献与阅读定位](#sec-11)

<a id="sec-1"></a>
## 1. 论文要解决什么问题

密码学中的困难实例通常由随机算法生成。“存在某些很难的格”本身并不能保证“随机生成的实例也很难”。这篇论文研究的，是如何把这两种困难性严格连接起来。

**严格地说，归约的方向是：**

$$
\text{最坏情况下的格问题}
\;\leq_{\mathrm{rand\text{-}poly}}\;
\text{均匀随机实例上的 SIS}.
$$

这里的含义是：**假设**有一个算法能以非忽略概率解决随机 SIS 实例，就能利用它构造随机多项式时间算法，解决相应的最坏情况格问题。反过来，在相应最坏情况困难性假设下，可以推出随机 SIS 的困难性。论文并没有无条件证明格问题难，也没有给出一个现实中可直接调用的高效 SIS 求解器。

本文使用欧氏范数版本的 SIS。给定

$$
A\in\mathbb Z_q^{\,n\times m},
$$

寻找

$$
z\in\mathbb Z^m\setminus\{0\},
\qquad
Az=0\pmod q,
\qquad
\|z\|_2\leq\beta.
$$

平均情况是指 $A$ 在整个 $\mathbb Z_q^{\,n\times m}$ 上均匀采样。参数 $m,q,\beta$ 必须满足相应定理的条件；不能把“随机矩阵”和“任意短向量界”拼在一起就引用困难性。

论文用高斯分析把若干最坏情况格问题与平均情况 SIS 之间的连接因子做到近线性量级 $\widetilde O(n)$。其中 $\widetilde O$ 隐去对数因子。这是概括，不表示每个问题都具有完全相同的精确参数，也不表示得到精确最短向量算法。第 8 节会给出具体区分。

对初学者而言，先抓住以下两个证明任务：

- **分布任务：** 从任意输入格出发，生成接近均匀的 SIS 输入。
- **几何任务：** 将 SIS 返回的短整数关系转换成原格中有用的短向量或解码结果。

平滑性质主要处理第一个任务；条件离散高斯及其矩估计是第二个任务的重要支撑。

<a id="sec-2"></a>
## 2. 符号与基础对象

### 2.1 格、基本区域和模格操作

全文令 $\Lambda=B\mathbb Z^n\subset\mathbb R^n$ 为满秩格，$B$ 的列向量构成一组基。记

$$
d=\det(\Lambda)=|\det B|,
\qquad
P(B)=B[0,1)^n.
$$

$P(B)$ 是基本平行多面体，体积为 $d$。每个 $x\in\mathbb R^n$ 都能唯一写成

$$
x=Bk+y,\qquad k\in\mathbb Z^n,\quad y\in P(B).
$$

具体地，

$$
x\bmod B=x-B\lfloor B^{-1}x\rfloor.
$$

向量的下取整逐坐标进行。这个操作选取商空间 $\mathbb R^n/\Lambda$ 的一个代表元；它不是“寻找距离 $x$ 最近的格点”。

**直观解释。** 可以把空间中所有基本区域平移叠到同一个 $P(B)$ 上。基改变时代表区域会改变，但被识别为同一个点的关系 $x-y\in\Lambda$ 不变。

### 2.2 连续高斯与离散高斯

采用论文的参数约定：

$$
\rho_{s,c}(x)
=\exp\!\left(-\pi\frac{\|x-c\|^2}{s^2}\right),
\qquad
\rho_s=\rho_{s,0},
\qquad s>0.
$$

对于可数集合 $S$，记 $\rho_{s,c}(S)=\sum_{x\in S}\rho_{s,c}(x)$。$\rho$ 是权重函数，不是已经归一化的概率密度。

连续高斯 $D_{s,c}$ 的密度是

$$
f_{s,c}(x)=\frac{\rho_{s,c}(x)}{s^n},
\qquad
\int_{\mathbb R^n}\rho_{s,c}(x)\,dx=s^n.
$$

其均值为 $c$，协方差为 $\frac{s^2}{2\pi}I_n$。因此每个坐标的标准差是 $s/\sqrt{2\pi}$，**不是 $s$**；离中心的均方根距离是 $s\sqrt{n/(2\pi)}$。

格上的离散高斯定义为

$$
\Pr_{X\sim D_{\Lambda,s,c}}[X=x]
=\frac{\rho_{s,c}(x)}{\rho_{s,c}(\Lambda)},
\qquad x\in\Lambda.
$$

它在格点上按高斯权重分配概率。将连续高斯样本直接舍入到最近格点，一般不会得到这个精确分布。

### 2.3 统计距离与可忽略误差

连续密度 $p,q$ 的统计距离，也称总变差距离，为

$$
\Delta(P,Q)=\frac12\int |p(x)-q(x)|\,dx.
$$

离散情形把积分换成求和。它控制任意事件的概率差，因此也控制任意算法在两种输入分布上的成功概率差。

如果各列独立且每列距理想分布至多 $\delta$，则 $m$ 列联合分布的统计距离至多 $m\delta$。对两种分布执行相同的确定性或随机处理，不会增大统计距离。

密码学中的“可忽略”是一个关于维数／安全参数 $n$ 的渐近概念：

$$
\varepsilon(n)\text{ 可忽略}
\quad\Longleftrightarrow\quad
\forall C>0,\ 
\varepsilon(n)<n^{-C}
\text{ 对充分大的 }n\text{ 成立}.
$$

固定很小的常数不是可忽略函数；固定指数的 $n^{-100}$ 也不是。

<a id="sec-3"></a>
## 3. 对偶格、Fourier 与 Poisson 各负责什么

### 3.1 对偶格给出与格周期相容的频率

对偶格定义为

$$
\Lambda^*
=\{w\in\mathbb R^n:\langle w,\lambda\rangle\in\mathbb Z
\text{ 对所有 }\lambda\in\Lambda\},
\qquad
\Lambda^*=B^{-T}\mathbb Z^n.
$$

考虑复指数波

$$
\chi_w(x)=e^{2\pi i\langle w,x\rangle}.
$$

它具有格周期，当且仅当

$$
\chi_w(x+\lambda)=\chi_w(x)
\quad\forall\lambda\in\Lambda
\quad\Longleftrightarrow\quad
w\in\Lambda^*.
$$

理由是两者的比值为 $e^{2\pi i\langle w,\lambda\rangle}$，而 $e^{2\pi it}=1$ 当且仅当 $t\in\mathbb Z$。

**直观解释。** 对偶格记录了哪些“波”能在基本区域的边界处接起来。这里有所有 $w\in\Lambda^*$ 对应的频率，并非只有 $n$ 个所谓“基本波”。$w=0$ 对应常数，其他频率描述周期密度的不均匀起伏。

例如，

$$
\Lambda=2\mathbb Z\times3\mathbb Z
\quad\Longrightarrow\quad
\Lambda^*=\tfrac12\mathbb Z\times\tfrac13\mathbb Z.
$$

原格在某方向的周期越长，相应的允许频率就越低。低频在高斯平滑后衰减较慢，这是短对偶向量与平滑参数联系紧密的原因。

### 3.2 Fourier 变换把高斯宽度变成频率衰减

采用

$$
\widehat f(w)=
\int_{\mathbb R^n}f(x)e^{-2\pi i\langle x,w\rangle}\,dx.
$$

对应的高斯变换为

$$
\widehat{\rho_{s,c}}(w)
=s^n e^{-\pi s^2\|w\|^2}e^{-2\pi i\langle c,w\rangle}.
$$

除以归一化因子 $s^n$ 后，

$$
\widehat f_{s,c}(w)
=e^{-\pi s^2\|w\|^2}e^{-2\pi i\langle c,w\rangle}.
$$

需要记住两点：

- 空间中的高斯越宽，非零频率的系数 $e^{-\pi s^2\|w\|^2}$ 越小。
- 改变中心 $c$ 只改变相位，不改变系数的绝对值。

### 3.3 Poisson 求和连接原格与对偶格

对 Schwartz 函数 $f$，Poisson 求和公式给出

$$
\sum_{\lambda\in\Lambda}f(\lambda)
=\frac1d\sum_{w\in\Lambda^*}\widehat f(w).
$$

Schwartz 条件保证函数及其各阶导数足够快地衰减。本文使用的高斯以及“多项式乘高斯”都满足条件，因而不需要对任意函数不加条件地套用公式。

其作用是把**原格上的空间求和**改写成**对偶格上的频率求和**。特别是对高斯，右侧各项具有显式的指数衰减，容易将零频率主项与非零频率误差分开。

这三个工具的分工是：对偶格确定求和频率，Fourier 变换计算每个频率的权重，Poisson 求和把它们与需要研究的格上分布连接起来。

<a id="sec-4"></a>
## 4. 平滑参数：定义、意义与精确结论

### 4.1 严格定义

**原文 Definition 3.1，p. 10。** 对 $\varepsilon>0$，

$$
\boxed{
\eta_\varepsilon(\Lambda)
=\min\left\{
s>0:
\rho_{1/s}(\Lambda^*\setminus\{0\})\leq\varepsilon
\right\}.
}
$$

展开后就是

$$
\rho_{1/s}(\Lambda^*\setminus\{0\})
=\sum_{w\in\Lambda^*\setminus\{0\}}
e^{-\pi s^2\|w\|^2}.
$$

这里的 $1/s$ 很关键：$\rho_t(w)=e^{-\pi\|w\|^2/t^2}$，代入 $t=1/s$ 才得到上式。$s$ 增大时，该和减小。

这个和是连续、严格递减的函数；在 $s\to0^+$ 时趋于无穷，在 $s\to\infty$ 时趋于零。所以定义中的最小值存在。参数 $\varepsilon$ 越小，要求越严格，$\eta_\varepsilon$ 越大。

**直观解释。** $\eta_\varepsilon$ 是一个高斯宽度门槛：超过它以后，与格周期相容的所有非零频率的**绝对权重总和**至多为 $\varepsilon$。控制单个频率很小不够，因为必须把所有频率累加起来。

“高斯参数 $s$”是我们选取的宽度；“平滑参数 $\eta_\varepsilon(\Lambda)$”是由格与误差要求共同决定的门槛。它们不是两个名称不同的同一个对象。

### 4.2 核心结论：连续高斯模格接近均匀

**严格结论（原文 Lemma 4.1，p. 12）。** 令 $X\sim D_{s,c}$，$Y=X\bmod B$。对任意 $s>0$ 与 $c\in\mathbb R^n$，

$$
\Delta\bigl(Y,U(P(B))\bigr)
\leq
\frac12\rho_{1/s}(\Lambda^*\setminus\{0\}).
$$

特别地，$s\geq\eta_\varepsilon(\Lambda)$ 时，

$$
\boxed{
\Delta\bigl(D_{s,c}\bmod B,U(P(B))\bigr)\leq\varepsilon/2.
}
$$

下面完整推导这个结论。

**第一步：写出折叠后的密度。** 对 $y\in P(B)$，所有 $y+\lambda$ 都映到同一个代表点，因此

$$
p_{s,c}(y)
=\sum_{\lambda\in\Lambda}f_{s,c}(y+\lambda)
=\frac1{s^n}\sum_{\lambda\in\Lambda}\rho_{s,c}(y+\lambda).
$$

**第二步：对平移高斯使用 Poisson 求和。** 得到

$$
p_{s,c}(y)
=\frac1d
\sum_{w\in\Lambda^*}
e^{-\pi s^2\|w\|^2}
e^{2\pi i\langle w,y-c\rangle}.
$$

**第三步：分离零频率。** 均匀分布的密度为 $1/d$，而 $w=0$ 恰好贡献 $1/d$：

$$
p_{s,c}(y)-\frac1d
=\frac1d
\sum_{w\in\Lambda^*\setminus\{0\}}
e^{-\pi s^2\|w\|^2}
e^{2\pi i\langle w,y-c\rangle}.
$$

由三角不等式，

$$
\left|p_{s,c}(y)-\frac1d\right|
\leq\frac1d\rho_{1/s}(\Lambda^*\setminus\{0\}).
$$

当 $s\geq\eta_\varepsilon(\Lambda)$ 时，得到逐点界

$$
\boxed{
\frac{1-\varepsilon}{d}
\leq p_{s,c}(y)\leq
\frac{1+\varepsilon}{d}.
}
$$

下界在 $0<\varepsilon<1$ 时尤其有用。

**第四步：积分。** 因 $\operatorname{vol}(P(B))=d$，

$$
\begin{aligned}
\Delta\bigl(Y,U(P(B))\bigr)
&=\frac12\int_{P(B)}
\left|p_{s,c}(y)-\frac1d\right|dy\\
&\leq\frac12\cdot d\cdot\frac{\varepsilon}{d}
=\frac{\varepsilon}{2}.
\end{aligned}
$$

这就是定义中对偶格高斯和的作用：它直接控制密度相对于均匀密度的起伏。

### 4.3 “平滑”究竟平滑了什么

**直观解释。** 在每个格点放置一个相同的高斯形状，再将这些形状叠加。宽度小时还能看到峰谷，宽度增大后峰谷差缩小。上面的周期密度公式就是这幅图像的严格版本。

但要避免三个误解：

1. **不是把无限空间变成均匀概率分布。** $\mathbb R^n$ 和无限格都不存在这样的均匀概率分布。严格结论位于有限体积基本区域，等价地位于商空间 $\mathbb R^n/\Lambda$。
2. **不是让离散高斯变成“全格均匀分布”。** 模格接近均匀的对象是连续高斯经模格后的分布。
3. **不是几乎处处完全一样。** 有限的 $s$ 下通常仍有非零频率；结论是误差可控的近似均匀。

**本文推导：它也是一个精确的逐点误差门槛。** 令

$$
R_s=\rho_{1/s}(\Lambda^*\setminus\{0\}).
$$

上面的估计不仅给出上界，而且

$$
\sup_{y\in P(B)}|d\,p_{s,c}(y)-1|=R_s.
$$

取 $y=c\bmod B$ 时，$y-c\in\Lambda$，所有相位都等于 $1$，故达到上界。因此 $s\geq\eta_\varepsilon$ 等价于上述逐点相对误差不超过 $\varepsilon$。不过，仅知道统计距离不超过 $\varepsilon/2$，不能倒推出同样的平滑参数条件：积分误差和最大误差是不同标准。

### 4.4 一个一维例子

取 $\Lambda=a\mathbb Z$，$a>0$，则 $\Lambda^*=a^{-1}\mathbb Z$。对中心 $c=0$，

$$
p_s(y)=\frac1a\left[
1+2\sum_{k=1}^{\infty}
e^{-\pi(s/a)^2k^2}
\cos\!\left(\frac{2\pi ky}{a}\right)
\right],
\quad 0\leq y<a.
$$

平滑参数控制的量是

$$
R_s=2\sum_{k=1}^{\infty}e^{-\pi(s/a)^2k^2}.
$$

令 $u=\pi(s/a)^2$，利用 $k^2\geq k$，得到

$$
2e^{-u}\leq R_s
\leq2\sum_{k=1}^{\infty}e^{-uk}
=\frac{2}{e^u-1}.
$$

因而在 $0<\varepsilon<1$ 时，

$$
a\sqrt{\frac{\ln(2/\varepsilon)}{\pi}}
\leq\eta_\varepsilon(a\mathbb Z)
\leq
a\sqrt{\frac{\ln(1+2/\varepsilon)}{\pi}}.
$$

这是**本文推导的严格上下界**，不是一般格上的精确公式。它说明一维中门槛随格间距线性变化，而对误差要求的依赖是对数量级的平方根。

<a id="sec-5"></a>
## 5. 怎样用格的几何量控制平滑参数

平滑参数定义在对偶格上。若不能把它与熟悉的格几何量联系起来，归约最后就难以表述成标准格问题。

记 $\lambda_i(\Lambda)$ 为第 $i$ 个逐次极小值，即最小的 $r$，使半径为 $r$ 的原点中心球内含有 $i$ 个线性无关格向量。特别地，$\lambda_1$ 是最短非零格向量长度，$\lambda_n$ 控制一组 $n$ 个独立格向量的最大长度。这组向量不一定是整格的一组基。

### 5.1 最短对偶向量揭示必要尺度

取最短非零对偶向量 $w$。$w$ 和 $-w$ 都贡献权重，所以

$$
R_s\geq2e^{-\pi s^2\lambda_1(\Lambda^*)^2}.
$$

于是，对于 $0<\varepsilon<1$，

$$
\eta_\varepsilon(\Lambda)
\geq
\frac{\sqrt{\ln(2/\varepsilon)/\pi}}
{\lambda_1(\Lambda^*)}.
$$

这是本文从定义推得的下界。它说明：很短的对偶向量对应不容易被压低的低频，要求更宽的高斯。

**原文 Lemma 3.2，p. 10** 还给出特定误差下的上界：

$$
\eta_{2^{-n}}(\Lambda)
\leq \frac{\sqrt n}{\lambda_1(\Lambda^*)}.
$$

其证明利用 Banaszczyk 的格高斯尾界：将对偶格放大，使所有非零格点位于半径 $\sqrt n$ 的球外，再控制这些点的总高斯权重。这里的误差是特定的 $2^{-n}$，不能直接把下标改成任意更小的 $\varepsilon$。

### 5.2 与 $\lambda_n$ 的联系

**严格结论（原文 Lemma 3.3，p. 11）。** 对任意 $\varepsilon>0$，

$$
\boxed{
\eta_\varepsilon(\Lambda)
\leq
\sqrt{\frac{\ln(2n(1+1/\varepsilon))}{\pi}}\,
\lambda_n(\Lambda).
}
$$

下面解释证明的关键结构。选取独立的 $v_1,\ldots,v_n\in\Lambda$，使 $\|v_i\|\leq\lambda_n(\Lambda)$，定义对偶格的切片

$$
S_{i,j}
=\{w\in\Lambda^*:\langle w,v_i\rangle=j\},
\qquad j\in\mathbb Z.
$$

由于内积是整数，每个切片离零超平面的距离为 $|j|/\|v_i\|$。对于非空切片，先将它沿法向移回零超平面，再使用“格高斯质量在格点中心最大”（原文 Lemma 2.9），得到

$$
\rho_{1/s}(S_{i,j})
\leq
e^{-\pi s^2j^2/\|v_i\|^2}
\rho_{1/s}(S_{i,0}).
$$

该比较在相应超平面内的格上进行；若切片为空，不等式显然成立。

令 $H=\rho_{1/s}(\Lambda^*)$、$T_i=\rho_{1/s}(\Lambda^*\setminus S_{i,0})$ 和 $a=\pi(s/\lambda_n)^2$。用几何级数估计，

$$
T_i\leq\frac{2}{e^a-1}(H-T_i),
\qquad
T_i\leq\frac{2H}{e^a+1}.
$$

因为 $v_1,\ldots,v_n$ 张成整个空间，任何非零 $w$ 都不可能同时与它们正交。因此对非零对偶格点作并集估计，记其总权重为 $R_s$，就有

$$
R_s\leq\sum_{i=1}^nT_i
\leq\frac{2n}{e^a+1}(1+R_s).
$$

取 $e^a=2n(1+1/\varepsilon)$，移项可得

$$
R_s\leq
\frac{2n}{e^a+1-2n}
<\varepsilon.
$$

这便给出了引理的宽度上界。

**直观解释。** 每个短原格向量都迫使对偶高斯质量主要集中到一个正交超平面；$n$ 个独立方向对应的超平面交集只有原点。于是非零对偶点的总权重很小。

### 5.3 如何读渐近记号

对指定的 $\varepsilon(n)$，应该优先使用上面的显式对数界。例如：

- $\varepsilon=n^{-C}$、$C>0$ 固定时，得到 $O(\sqrt{\log n})\lambda_n$ 的上界，但误差仍不是可忽略函数。
- $\varepsilon=e^{-(\ln n)^2}$ 时，误差可忽略，上界为 $O(\log n)\lambda_n$。
- 对任意 $h(n)=\omega(\sqrt{\log n})$，可选择某个可忽略的 $\varepsilon(n)$，使 $\eta_{\varepsilon(n)}(\Lambda)\leq h(n)\lambda_n(\Lambda)$ 对所有 $n$ 维格成立。

最后一句的量词是“给定这样的 $h$，存在合适的可忽略误差”，不是“对每一种可忽略误差都适用同一个 $h$”。例如把误差指定为 $2^{-n}$，上述 $\lambda_n$ 界一般只直接给出 $O(\sqrt n)\lambda_n$。

另外，由对偶格的缩放关系可得

$$
\eta_\varepsilon(a\Lambda)=a\,\eta_\varepsilon(\Lambda),
\qquad a>0.
$$

因此，平滑参数具有长度量纲，并且与基的选择无关。

<a id="sec-6"></a>
## 6. 为什么还需要离散高斯的矩估计

模格接近均匀，只说明生成的随机实例具有正确分布。它还没有说明将多个样本按 SIS 的解组合后，输出向量会足够短。

### 6.1 归一化因子稳定不等于分布相同

Poisson 求和还给出

$$
\rho_{s,c}(\Lambda)
=\frac{s^n}{d}
\sum_{w\in\Lambda^*}
e^{-\pi s^2\|w\|^2}e^{-2\pi i\langle w,c\rangle}.
$$

所以当 $0<\varepsilon<1$ 且 $s\geq\eta_\varepsilon(\Lambda)$ 时，

$$
(1-\varepsilon)\frac{s^n}{d}
\leq\rho_{s,c}(\Lambda)
\leq(1+\varepsilon)\frac{s^n}{d}.
$$

这个结论对任意中心 $c$ 都成立，说明离散高斯的归一化因子随中心改变时变化很小。

不过，不能说“离散高斯与连续高斯的统计距离很小”。把二者视为 $\mathbb R^n$ 上的概率测度，离散高斯以概率 $1$ 落在 $\Lambda$，连续高斯以概率 $0$ 落在 $\Lambda$，因此总变差距离为 $1$。后面的结论比较的是**某些矩**，而不是原始分布的总变差距离。

### 6.2 原文的一阶与二阶矩结论

**严格结论（原文 Lemma 4.2，pp. 12–14）。** 设

$$
X\sim D_{\Lambda,s,c},\quad
0<\varepsilon<1,\quad
s\geq2\eta_\varepsilon(\Lambda),\quad
\|u\|=1.
$$

则

$$
\left|\mathbb E\langle X-c,u\rangle\right|
\leq\frac{\varepsilon s}{1-\varepsilon},
$$

$$
\left|
\mathbb E\!\left[\langle X-c,u\rangle^2\right]
-\frac{s^2}{2\pi}
\right|
\leq\frac{\varepsilon s^2}{1-\varepsilon}.
$$

这是围绕指定中心 $c$ 的矩，第二个式子不是直接写出的中心化方差，因为 $\mathbb EX$ 不一定恰好等于 $c$。

**注意假设：这里是 $s\geq2\eta_\varepsilon$。** Lemma 4.1 的近均匀结论只要求 $s\geq\eta_\varepsilon$，不能不加说明地将两条引理的条件混用。

沿正交坐标方向求和，可得原文 Lemma 4.3：

$$
\|\mathbb E(X-c)\|^2
\leq
\left(\frac{\varepsilon}{1-\varepsilon}\right)^2s^2n,
$$

$$
\mathbb E\|X-c\|^2
\leq
\left(\frac1{2\pi}+\frac{\varepsilon}{1-\varepsilon}\right)s^2n.
$$

这说明位移的均值很小，其均方长度仍接近连续高斯的尺度。

### 6.3 证明思路：零频率给出连续高斯的矩

先将格与中心除以 $s$，归一化为 $s=1$；再旋转坐标，使 $u=e_1$。设

$$
g_1(x)=(x_1-c_1)\rho_{1,c}(x),
\qquad
g_2(x)=(x_1-c_1)^2\rho_{1,c}(x).
$$

所求矩是

$$
\frac{\sum_{\lambda\in\Lambda}g_j(\lambda)}
{\rho_{1,c}(\Lambda)}.
$$

对分子、分母使用 Poisson 求和。通过高斯导数与 Fourier 变换的关系，得到

$$
\widehat g_1(w)=-iw_1\widehat{\rho_{1,c}}(w),
$$

$$
\widehat g_2(w)=
\left(\frac1{2\pi}-w_1^2\right)
\widehat{\rho_{1,c}}(w).
$$

在 $w=0$ 处，这两个表达式恰好给出连续高斯的矩 $0$ 和 $1/(2\pi)$。剩下的任务是控制非零频率。

分母的归一化频率和至少为 $1-\varepsilon$。分子误差中出现

$$
|w_1|e^{-\pi\|w\|^2}
\quad\text{或}\quad
w_1^2e^{-\pi\|w\|^2}.
$$

这些多项式因子可以被较慢衰减的高斯 $e^{-\pi\|w\|^2/4}$ 吸收。归一化后的条件 $\eta_\varepsilon(\Lambda)\leq1/2$ 保证这个较宽的对偶高斯在非零点上的总权重仍至多为 $\varepsilon$。恢复尺度后，一阶矩乘 $s$，二阶矩乘 $s^2$，得到上述界。

这说明为什么“矩稳定性”比“归一化因子稳定性”需要额外论证：计算矩会引入频率的多项式因子。

<a id="sec-7"></a>
## 7. 高斯如何生成归约需要的样本

### 7.1 一个同时保留近均匀性与格点关联的采样方法

**严格结论（原文 Lemma 5.7，p. 20；理想实数采样描述）。** 给定格基 $B$、目标偏移 $t$ 和 $s\geq\eta_\varepsilon(\Lambda)$：

1. 采样连续高斯 $R\sim D_{s,t}$。
2. 令 $C=(-R)\bmod B$。
3. 令 $Y=R+C$。

于是 $C\in P(B)$、$Y\in\Lambda$，并且

$$
\Delta(C,U(P(B)))\leq\varepsilon/2.
$$

更关键的是，固定 $C=c$ 后，

$$
Y\mid C=c\ \sim\ D_{\Lambda,s,t+c}.
$$

**为什么条件分布是离散高斯？** 对每个 $y\in\Lambda$，与 $(C,Y)=(c,y)$ 相应的连续噪声是 $R=y-c$，其密度正比于

$$
e^{-\pi\|y-c-t\|^2/s^2}.
$$

对所有 $y\in\Lambda$ 的权重归一化，恰好得到以 $t+c$ 为中心的离散高斯。

这里 $C=c$ 是零概率事件，不能把普通条件概率的分子、分母直接相除。严格做法是使用连续坐标 $c$ 与离散坐标 $y$ 的联合密度，得到一个正则条件分布；上述公式提供了它的一个版本。算法也不是不断采样直到连续点恰好落在格上。

### 7.2 “能输出条件分布”不等于“能在任意给定中心采样”

这个程序输出的是一对相关变量 $(C,Y)$。它没有让我们预先任意指定 $C=c$ 后，再高效采样对应的 $Y$。

因此，Lemma 5.7 并不是 GPV 的短基离散高斯采样算法，也不是一个仅输入任意基、任意目标中心，就能在平滑参数尺度上独立采样的通用黑箱。这一点能避免把两篇论文的贡献混在一起。

<a id="sec-8"></a>
## 8. 从这些工具到最坏情况—平均情况归约

本节先给出便于检查的简化代数模型，再说明真正归约补上了什么。**简化模型不等于完整的格问题求解算法。**

### 8.1 近均匀的基本区域样本变成随机模方程

暂时以 $P(B)$ 为分区区域，将其分成 $q^n$ 个等体积小平行多面体。对独立样本 $c_1,\ldots,c_m$，定义

$$
a_i=\lfloor qB^{-1}c_i\rfloor\in\{0,\ldots,q-1\}^n,
\qquad
\widetilde c_i=\frac{B a_i}{q}.
$$

如果 $c_i$ 精确均匀，则标签 $a_i$ 在 $\mathbb Z_q^n$ 上精确均匀。若每个 $c_i$ 与均匀分布的统计距离至多为 $\varepsilon/2$，则

$$
\Delta\!\left(A,U(\mathbb Z_q^{\,n\times m})\right)
\leq m\varepsilon/2,
\qquad A=[a_1,\ldots,a_m].
$$

因此，若假设 oracle 在真正均匀输入上成功率为 $\delta(n)$，它在这些输入上的成功率至少为 $\delta(n)-m\varepsilon/2$。选取可忽略误差、且 $m$ 为多项式，才能保留非忽略成功概率。

### 8.2 SIS 解为什么能回到原格

oracle 返回 $z\neq0$ 且 $Az=0\bmod q$。于是 $Az/q\in\mathbb Z^n$，从而

$$
\widetilde C z=\frac{B Az}{q}\in\Lambda,
\qquad
\widetilde C=[\widetilde c_1,\ldots,\widetilde c_m].
$$

另一方面，每个 $y_i\in\Lambda$，所以 $Yz\in\Lambda$。因此

$$
v=(\widetilde C-Y)z\in\Lambda.
$$

它的长度可分成两部分：

$$
v=(\widetilde C-C)z+(C-Y)z.
$$

第一项来自将连续点映到小区域角点时的取整误差；第二项来自高斯位移的线性组合。

若 $\|B\|_{\max}=\max_i\|b_i\|_2$，则

$$
\|\widetilde c_i-c_i\|
\leq\frac{n\|B\|_{\max}}q,
$$

因此

$$
\|(\widetilde C-C)z\|
\leq\frac{n\|B\|_{\max}}q\|z\|_1
\leq\frac{n\sqrt m\,\beta}{q}\|B\|_{\max}.
$$

实际归约使用一组较短的独立格向量 $S$，并通过 Lemma 5.8 将样本适当地转移到 $P(S)$，使这里的误差取决于 $\|S\|_{\max}$。$S$ 不一定生成整个 $\Lambda$，所以这次转移不是无条件替换一个符号；需要原文的组合程序。

### 8.3 为什么不能直接把组合当作独立高斯之和

若 $z$ 预先固定，且 $e_i$ 是独立、均值为零的向量，则

$$
\mathbb E\left\|\sum_i z_i e_i\right\|^2
=\sum_i z_i^2\mathbb E\|e_i\|^2.
$$

这让长度尺度依赖于 $\|z\|_2$，而不是粗略三角不等式中的 $\|z\|_1$。

但在归约中，**$z$ 是 oracle 看过 $A$ 后选出的**，而 $A$ 又来自样本。不能把 $z$ 与所有原始随机量当成独立变量。

正确处理方法是：固定完整的 $C$、组合程序和 oracle 所见的信息及其输出 $z$，再分析没有公开给 oracle 的 $Y$。条件化后，各个 $y_i$ 仍按相应中心的离散高斯独立分布。于是可以使用第 6 节对任意中心均成立的矩界。

更一般地，若独立的 $e_i$ 满足

$$
\mathbb E\|e_i\|^2\leq L,
\qquad
\|\mathbb Ee_i\|^2\leq M,
$$

则对固定 $z$（原文 Lemma 2.11），

$$
\mathbb E\left\|\sum_i z_i e_i\right\|^2
\leq(L+mM)\|z\|_2^2.
$$

离散高斯的矩界给出 $L=O(s^2n)$、$M$ 为相对于 $s^2n$ 的可忽略量。因此在多项式数量的样本下，组合向量的均方长度为

$$
O(s^2n\beta^2).
$$

**这里是二阶矩估计，不是每次输出都满足的确定性长度界。** 进一步使用 Markov 不等式与成功概率分析，才得到归约需要的长度保证。

### 8.4 完整证明的中间任务：IncGDD

上面只说明“输出属于格且可能较短”。它尚未保证输出非零，更不能直接保证得到 $n$ 个独立向量。即使 $z\neq0$，也不能推出 $(\widetilde C-Y)z\neq0$。

原文引入增量保证距离解码问题 $\mathrm{IncGDD}^{\phi}_{\gamma,g}$。输入包括：

- 格基 $B$；
- $\Lambda$ 中 $n$ 个独立向量组成的 $S$；
- 目标点 $t$；
- 满足 $r>\gamma(n)\phi(\Lambda)$ 的数 $r$。

要求输出 $v\in\Lambda$，使

$$
\|v-t\|\leq\frac{\|S\|_{\max}}{g(n)}+r.
$$

这是相对于指定几何尺度 $\phi(\Lambda)$ 的距离保证，不是通常 CVP 中“相对最近距离的乘法近似”。

**原文 Theorem 5.9，pp. 21–24。** 在参数具有多项式位长、$m,\beta$ 多项式有界、$\varepsilon(n)$ 可忽略，且

$$
q(n)\geq g(n)n\sqrt{m(n)}\,\beta(n)
$$

时，可以将 $\phi=\eta_\varepsilon$、$\gamma(n)=\beta(n)\sqrt n$ 的 IncGDD 归约到平均情况 $\mathrm{SIS}_{q,m,\beta}$。

证明思路中的关键补全如下。

1. 随机猜测 oracle 解 $z$ 的一个非零坐标：猜 $j$ 和 $\alpha=z_j$。设 $t_j=-t/\alpha$，其余 $t_i=0$，则猜中时 $Tz=-t$。
2. 取高斯宽度 $s=2r/\gamma>2\eta_\varepsilon(\Lambda)$，使采样引理和矩界都适用。
3. 用组合程序得到 $x\in\Lambda$，满足 $\|x-Cz\|\leq\|S\|_{\max}/g$。输出 $v=x-Yz$。
4. 猜中且 oracle 成功时，
   $$
   \|v-t\|
   \leq\frac{\|S\|_{\max}}g+
   \|(Y-C-T)z\|.
   $$
   条件二阶矩分析给出第二项不超过 $r$ 的常数成功概率。猜坐标及 oracle 的成功概率只造成多项式损失，分布替换造成可忽略损失。
5. 后续最坏情况之间的归约再把 IncGDD 用于独立短向量、保证距离解码及覆盖半径问题。这些步骤承担了简化模型没有解决的任务。

上述过程也说明：平滑参数不用作为一个可以精确高效计算的量直接输入算法；在该中间问题中，输入的 $r$ 与承诺 $r>\gamma\eta_\varepsilon$ 已为宽度选择提供了依据。

### 8.5 近线性连接因子从哪里来

针对 SIVP 等问题，可以这样追踪量级：

$$
\text{组合噪声尺度}
\ \sim\ \beta\sqrt n\,\eta_\varepsilon(\Lambda).
$$

对于论文选择的 $m=\Theta(n\log n)$ 与合适多项式模数，可以取保证 SIS 解存在的

$$
\beta=\sqrt m\,q^{n/m}=O(\sqrt{n\log n}).
$$

再把第 5 节的平滑参数上界代入，得到近线性的连接尺度。选择例如 $\varepsilon=e^{-(\ln n)^2}$，这条链直接给出 $O(n(\log n)^{3/2})\lambda_n$ 量级；若仔细选择可忽略误差，可以得到作者稿中更精确的定理表述。

具体来说，**Theorem 5.16（p. 25）** 表明：对 $m(n)=\Theta(n\log n)$，存在合适的 $q(n)=O(n^2\log n)$，使任意 $\gamma(n)=\omega(n\log n)$ 的下列问题可归约到平均情况 SIS：

- $\mathrm{SIVP}_\gamma$：找 $n$ 个独立格向量，每个长度至多为 $\gamma\lambda_n$；
- $\mathrm{GDD}^{\lambda_n}_\gamma$：给定任意目标，找距其至多 $\gamma\lambda_n$ 的格点；
- $\mathrm{GapCRP}_\gamma$：覆盖半径的相应间隙判定问题。

这里的“存在合适的 $q$”不能改写成“每个满足大 O 上界的 $q$ 都可以”。

**GapSVP 是另一条证明分支。** 作者稿 §5.4 经由特殊的 $\mathrm{GapCVP}'$ 与 $\mathrm{SIS}'$，使用对偶格中的见证及额外的尾界、Fourier 估计。$\mathrm{SIS}'$ 要求解至少有一个奇数坐标；当 $q$ 为奇数时，可以从普通 SIS 的解中不断约去公共因子 $2$，得到它的解。

**Theorem 5.23（p. 28）** 的特化结果使用 $m=\Theta(n\log n)$、合适的 $q=O(n^{5/2}\log n)$，得到 $\gamma=O(n\sqrt{\log n})$ 的 GapSVP 连接因子。选取满足条件的奇数模数后可转为普通 SIS。GapSVP 是区分“$\lambda_1\leq d$”与“$\lambda_1>\gamma d$”的承诺判定问题，不能据此声称得到同因子的 search-SVP 归约。

> **未展开：** 本节没有重建 §5.4 的完整验证器、抽样次数及成功概率证明。因此不能把第 8.2 节那条简单的短向量组合链当成 GapSVP 的完整证明。此处的作用是交代定理范围和工具分工。

<a id="sec-9"></a>
## 9. 论文贡献与容易混淆的说法

### 9.1 主要贡献

**第一，建立可用于归约的平滑参数分析。** 非零对偶高斯质量同时连接了格的几何结构与模格分布的近均匀性。Definition 3.1、Lemma 3.3 与 Lemma 4.1 是理解这条联系的核心。

**第二，证明条件化后仍能利用高斯的矩性质。** 归约不能只证明标签近均匀，还必须在 oracle 已看过标签并选择整数关系后分析隐藏位移。离散高斯的一、二阶矩控制使这种分析成立。

**第三，改善最坏情况到平均情况的连接因子。** 论文将若干问题的连接因子推进到近线性量级，并以高斯采样取代此前较复杂的区域计数分析。精确结论需按问题分别读取，不能仅凭摘要中的 $\widetilde O(n)$ 抹去条件。

原文 §1.1 对较早方法的描述是：既希望分区足够小以控制位移，又希望每个分区含有相近数量的格点以保证标签近均匀。高斯方法利用模格后的密度直接控制等体积分区的概率，减轻了这两种要求之间的冲突。**这篇论文的归约仍包含将连续样本变成离散标签的取整步骤**，不能将它描述成“完全避免 rounding”。

### 9.2 几个应该保留的逻辑边界

- **近均匀不等于完整归约。** 还需要代数映射、长度控制、目标或非零性处理，以及成功概率和运行时间分析。
- **矩接近不等于分布接近。** 连续与离散高斯比较的是特定统计量，不能直接声称它们在 $\mathbb R^n$ 上总变差很小。
- **短对偶向量重要，但不能只检查一个频率。** 平滑参数控制的是所有非零频率权重之和。
- **$\lambda_n$ 不是“给定基中最长向量长度”。** 一个很差的基可以很长，而格本身的 $\lambda_n$ 不变。
- **小误差不是可忽略误差。** 最终是否可用，取决于误差、样本数与 oracle 成功率之间的量化关系。
- **oracle 是条件假设。** 最坏到平均归约的意义是传递困难性，并非证明某个具体实例无条件不可解。
- **本文不是 GPV。** 本文不把 SampleD、短基陷门、签名或 IBE 当作 Micciancio–Regev 这篇论文的结果。

把各模块连起来看，本文的理解是：格周期结构决定对偶频率；高斯使这些频率显式衰减；Poisson 求和把衰减转换为近均匀性和矩控制；这两种控制分别保证归约的随机输入足够像真正随机输入，以及返回的整数关系可以转化成有用的格几何信息。

<a id="sec-10"></a>
## 10. 核对记录与待进一步核实的部分

### 10.1 本稿已完成的核对

以下结论已对照 [2005-12-14 作者稿][paper]，并在本文中保持相同的参数约定：

- Definition 3.1：非零对偶高斯质量使用参数 $1/s$。
- Lemma 3.2：$\varepsilon=2^{-n}$ 时的 $\sqrt n/\lambda_1(\Lambda^*)$ 上界。
- Lemma 3.3：$\ln(2n(1+1/\varepsilon))/\pi$ 的平方根及 $\lambda_n$ 因子。
- Lemma 4.1：任意中心、逐点密度估计，以及采用 $\frac12 L^1$ 约定后的 $\varepsilon/2$ 统计距离。
- Lemmas 4.2–4.3：要求 $s\geq2\eta_\varepsilon$；二阶矩以 $s^2/(2\pi)$ 为连续高斯参照。
- Lemma 5.7：$C=(-R)\bmod B$ 与 $Y=R+C$ 的符号，以及条件中心 $t+c$。
- Theorem 5.9：$q\geq gn\sqrt m\,\beta$、$\gamma=\beta\sqrt n$，并区分增量解码与直接求非零短向量。
- Theorems 5.16、5.23：两条分支的模数、近似因子、问题版本与奇数模数条件。
- 本文附加的一维界、缩放关系与逐点误差等式已从所列公式推导；它们没有被冒充成原文的新定理。

### 10.2 待核实／未展开

这些条目标记的是本文的边界，不是否定上述已核对的结论：

- [ ] **期刊版本对照。** 尚未逐页比较 2007 期刊排版版与本次使用的 2005 作者稿。正式引用期刊版时，应核实编号和是否存在修改；当前笔记的定位一律以链接的作者稿为准。
- [ ] **Lemma 5.8 的完整实现。** 本文只使用组合程序的输入输出性质，没有逐行复核从 $P(B)$ 到子格基本区域 $P(S)$ 的随机映射及全部算法细节。若要将第 8 节扩展为完整归约讲义，需要补上这部分。
- [ ] **§5.4 的完整 GapSVP 证明。** 定理表述已核对；验证器的谱条件、特征函数估计、Hoeffding 界和所有失败事件的合并尚未在本文逐项重证。
- [ ] **有限精度实现。** 第 7–8 节按论文的理想连续采样模型讲解。若要实现或给出按位复杂度的完整算法证明，应另行核对高斯采样近似、实数表示、取整边界误差，以及它们累积后的统计距离。
- [ ] **历史与后续改进。** 本文只讨论这篇论文及其内部结论，没有系统比较此后所有 SIS 归约。文中的近线性结果不应被理解为“截至 2026 年最优参数”。

<a id="sec-11"></a>
## 11. 参考文献与阅读定位

1. Daniele Micciancio and Oded Regev. **Worst-case to Average-case Reductions based on Gaussian Measures.** [作者公开稿，2005-12-14，PDF][paper]。本文的主要数学依据与编号来源。
2. [作者论文主页][homepage]。期刊信息：*SIAM Journal on Computing*, 37(1):267–302, 2007；初版 FOCS 2004。
3. [DOI: 10.1137/S0097539705447360][doi]。期刊发表记录；本文没有声称已逐页核对该排版版本。

建议按以下顺序对照原文：

- **§2，pp. 5–9：** 分布、格与 Gaussian/Fourier/Poisson 的符号约定。
- **§3，pp. 10–11：** 平滑参数定义及与几何量的联系。
- **§4，pp. 12–15：** 先读 Lemma 4.1，再读矩估计；理解近均匀与条件矩分别解决什么问题。
- **§5.2，pp. 18–24：** 先读 Lemma 5.7 和 Theorem 5.9，再补 Lemma 5.8。
- **§5.3，pp. 24–26：** 看中间问题如何转回标准格问题。
- **§5.4，pp. 27–31：** 单独阅读 GapSVP 分支，不将它与 SIVP 分支混为一条证明。

> 更正说明：若发现公式、假设或逻辑遗漏，欢迎指出具体位置并附原文依据。后续修改应优先修正数学结论及其适用条件，再改进直观解释。

[paper]: https://cims.nyu.edu/~regev/papers/average.pdf
[homepage]: https://cseweb.ucsd.edu/~daniele/papers/Gaussian.html
[doi]: https://doi.org/10.1137/S0097539705447360



