# 模型、来源与限制

## 数学模型

样品 i 的精确相对年距为 d_i，起点为 T，单年参考均值为 C(T+d_i)。测年值 y_i 和共同偏移 delta 满足：

    y_i | T,delta ~ Normal(C(T+d_i)+delta, s_i²+c_i(T+d_i)²)
    delta ~ Normal(mu,tau²)

位置先验在用户显式搜索的整数年份内均匀。各测量误差独立，各参考年份误差也独立；不同待测点不能重用同一个参考年份。

对共同偏移解析积分后：

    y | T ~ MVNormal(C_T + mu*1, diag(s_i²+c_i²) + tau²*1*1')

内核用等价的对角项加秩一项计算，不抽样、不随机近似。保留正态密度的行列式/归一化项，避免不同年份参考误差变化造成偏差。tau=0 表示固定偏移。

偏移输出的总体均值与标准差由各年份条件分布的混合计算，包含年代不确定度。混合分布可能多峰，不能据此假设偏移后验为单一正态。

## 参考准备的假设

- IntChron 日历换算为 astronomical_year=1950-calage。只导入北半球树轮、样品年宽为 1、三项日历误差为 0 的行。
- 带非零 r_date_sigma_batch 的记录包含批次共有误差。首版没有批次协方差，因此排除并记录。额外误差、非单位误差倍率、储库修正亦不支持，不能忽略。
- Brehm 2021 PANGAEA 年龄/误差的 ka 换为年，日历 ka BP 换为整数年份；EL 按原发布说明解释为同一年早+晚材。重复靶以原标识保留。来源没有单独编码的相关误差仍可能存在，不据此声称批次独立性已获证明。
- 单年均值误差使用显式 Birge 扩大约定。没有推断地区偏移随时间的变化，也没有完整模拟参考记录之间的相关性。

## 文献与实现身份

- Bronk Ramsey, van der Plicht & Weninger (2001), ‘Wiggle Matching’ Radiocarbon Dates, Radiocarbon 43, 381–389：[作者机构全文](https://pure.rug.nl/ws/portalfiles/portal/3102767/2001RadiocarbonBronkRamsey.pdf)。本工具采用已知年距联合似然的思想，不声称复现其全部示例。
- Bronk Ramsey (2009), Dealing with Outliers and Offsets in Radiocarbon Dating：[原始论文](https://www.cambridge.org/core/journals/radiocarbon/article/dealing-with-outliers-and-offsets-in-radiocarbon-dating/500255B3C1ABE8B3FBE4614433C0D1B4)。首版只处理一个共同正态偏移，没有自动异常值模型。
- [OxCal 已知年距模型文档](https://c14.arch.ox.ac.uk/oxcalhelp/hlp_analysis_oper.html)；[历年定义](https://c14.arch.ox.ac.uk/oxcalhelp/hlp_analysis_inform.html)；[IntChron 元数据定义](https://intchron.org/schema.html)。外部核查日期 2026-11-01。
- Pearson 等 2026 PNAS 正文第 2–3、6 页为应用背景，N(0,50²) BP 是其设置之一。0.2公开运行接口移除默认偏移标准差，要求按样品与参考来源显式选择。原代码、完整参考数据与模型未复现；不要把本工具的输出描述为作者结果。
- [Bayliss et al. 2024](https://doi.org/10.1017/RDC.2024.54)提供罗马时期反例；本工具已运行43组相关窗口并观察到多项未覆盖，见validation.md。共同高斯偏移不能自动解决参考曲线在某年代的形状或时变偏差。

## 不能直接推断的内容

末轮缺失、边材估计、砍伐/建造/毁弃间隔、储库效应、年距不确定、跨半球修正、参考间协方差、时变偏移均不在首版模型内。没有可信的参考覆盖和样品年轮位置，就没有可解释的 wiggle-match。单年取样不保证单年定年。
