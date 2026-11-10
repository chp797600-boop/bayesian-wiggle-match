# 输入与运行

解释运行结果时按[五关判读](five_gate_diagnostics.md)检查计算有效性、统计一致性、唯一性与精度、偏移合理性及独立真值。0.4.0自动计算适用固定偏移的Acomb，未知共享偏移不套该公式；外部OxCal实跑仍是独立对照。HPD区间数不等于年代簇数，偏移收紧不自动证明物理解释成立。

唯一性报告联合检查主峰及竞争年代带、完整HPD结构，以及局部年窗概率。0.3.0自动生成连续1/2/5/10/20/50年窗口、候选局部峰/平顶、谷部和盆地质量，保留全部数值并列最优起点；不自动输出“唯一性通过”。局部峰包括尾部细小起伏，盆地是确定性分界的记账区域，不是已验证的年代簇。详见五关判读的第三关。

## 年代及单位

整数天文年：公元 n 年写 n；公元前 n 年写 1-n。因此 -1600 是 1601 BCE，不是 1600 BCE。年份起点对应 ring_offset=0 的最早测年轮，不自动对应整棵树的末轮或砍伐年。

年龄采用常规 14C BP，误差为 1 sigma。首版不接受 Delta14C、F14C 或 ka BP 直接输入年龄列；必须先进行有来源、保留原值及转换记录的单位整理。

## 样品 CSV

必需字段：

|字段|含义|
|---|---|
|sample_id|样品标签|
|ring_offset|距最早测年轮的精确年数；非负整数，最小为 0；不能重复|
|c14_age_bp|测得的常规 14C 年龄|
|sigma_bp|年龄的正值 1 sigma 误差|
|sample_span_years|固定为 1；不能把合测的多轮写为 1|
|measurement_ids|原测量标识，多个用分号分隔|
|tree_id|树木/保守来源组标识，不能与参考重叠|

同轮重复测量需先核对批次、独立性和合并条件，不能把它们作为多个独立 ring_offset。外部新数据应填真实实验室号/测量号及树木标识，不要为了通过独立性检查改名。

## 年度参考 CSV

每个 astro_year 只能一行，必需字段为 astro_year、c14_age_bp、sigma_bp、member_ids、tree_ids；标识列表用分号分隔。参考均值误差必须包含建立参考时决定保留的真实不确定度。相邻年度误差相关性、地区/季节差异仍需另行评估，首版不建模它们。

## 调用

在技能目录中，以当前可用 Python 执行：

```text
python -B scripts/run.py --samples assets/example_samples.csv --reference assets/example_reference.csv --start-min 350 --start-max 431 --offset-mean 0 --offset-sd 50 --out <项目内新建结果目录>
```

示例来自独立树木已知年代资料，样品文件有 30 个年度点，真起点为 400 CE。示例是多峰后验演示：MAP 并不等于真实年份。以实际文件和验证记录为准，不将示例结果当作新样品答案。

`--out` 必须尚不存在。命令会先核对输入和计算，再创建结果目录，不覆盖已有结果。`--allow-incomplete` 仅在知情的探索分析中显式添加；输出条件于可评估位置，不给参考缺口分配零年代概率。

0.2要求显式`--offset-sd`，省略时拒绝运行；Python的execute接口亦如此。上例50 BP只复现示例参数，不能直接作为新样品的先验依据。

输出：posterior.csv 保存每年概率、该年份条件下的共同偏移均值/标准差和对数边际似然；fit_at_map.csv 保存 MAP 对齐及条件偏移后的残差；summary.json 保存 HPD、偏移总体矩及警告；manifest.json 保存输入 SHA-256、程序版本/哈希及参数。

0.3.0额外输出diagnostics.json、candidate_peaks.csv、peak_valleys.csv；summary.json也包含posterior_diagnostics，manifest包含诊断模块散列。没有峰间谷时，peak_valleys.csv仅有表头。sensitivity各offset_sd子目录同样保留这些输出。

0.4.0另保存model_checks.json（statistical_consistency与offset_joint），summary也包含两项，manifest记录model_checks代码散列。固定偏移Acomb的参考归一化域和截断检查、偏移等尾区间、收紧倍数及关联须一并解释，不只报一个通过值。逐轮留一/PPC/搜索参考比较/模拟等较重计算由[assessment入口](assessment.md)显式执行，不在每个常规run里自动重跑。

只读分析旧的本地年度后验（不重新拟合）：

```text
python -B scripts/diagnostics.py --posterior <旧结果目录>/posterior.csv --summary <旧结果目录>/summary.json --out <项目内新建诊断目录>
```

输入必须是本地整数逐年概率质量及对应summary；程序核查归一化、搜索范围、覆盖和缺年列表。不是直接读取OxCal任意分辨率的JS密度。输出目录不得已存在，原输入不改写。

windows中的conditional_best_mass是已评估位置条件后验的最优连续窗口质量；maximizing_start_ranges为所有并列最优窗口起点范围，长度w的终点为起点+w−1。并列容差1e-12是数值比较容差，不是科学等效阈值。fully_assessed_best_mass另给窗口内年份全部可评估的最优值，仍条件于同一个已评估后验；不能因此宣布整个搜索范围已覆盖。窗口大于搜索范围或不存在完整窗口时相应值为null，不向窗口外补概率。候选盆地按谷部最低平台的中点向下取整分界，分界年归左侧，质量不可当作经验证的模式概率。

HPD 按离散年份概率排序构造，保留不连续区间及实际覆盖概率。等概率年份同时保留，实际质量可能大于标称水平。不能把多个区间的最小值至最大值改写成一个连续的 95% 区间。

## 多先验对照与已知年代评分

```text
python -B scripts/sensitivity.py --samples assets/example_samples.csv --reference assets/example_reference.csv --start-min 350 --start-max 431 --offset-sds 20 50 100 --out <项目内新建对照目录>
```

必须显式给出至少两个不同的非负偏移标准差。同一输入和搜索窗运行所有先验，输出各自完整结果及sensitivity.csv/json，不选择“最佳先验”。输入散列随报告保存。参考缺年仍默认拒绝；确需探索时显式添加`--allow-incomplete`。

报告中的hpd95_year_count为所有HPD区间内整数年份数之和，hpd95_interval_count为不连续区间数；map_span_years仅是不同先验的MAP跨度，不是可信区间。总变差距离为0.5×sum(abs(p_a(t)-p_b(t)))，范围0—1，在相同可评估位置计算，表示先验改变导致的分布差异，不是准确率或科学适用性分数。即使所有先验MAP相同，整个分布也可能明显不同。

仅检验已有独立年代时，可加`--truth-start 400`。真值只在模型拟合后评分，不传给匹配函数、不改变后验。truth_status区分covered、missed、unassessed、outside_search_window及not_supplied；后三种truth_in_95_hpd为null，不应被当成False加入覆盖失败数。生产未知样品不要填写猜测的真值。

## 同轮合并与参考血缘

prepare.build_reference的Birge扩大是本工具既有参考聚合约定；它不等于所有论文对同轮重复测量采用的普通逆方差标准误。复现已发表案例时先核实原文的合并规则，保留两次原始测量、合并公式及输出。不要为匹配期望结果改动误差；必要的口径对照必须单独列出，不能覆盖主检验。

同一条已知年表可能跨多棵树。本轮罗马检验的跨树精确年距已由原文年表提供，不代表未知浮动木材也有这些信息。对于模型校准曲线，tree_ids不能代替底层输入数据的血缘审查；已知有同树输入时即使测量批次不同，也必须披露。

## 整理已下载资料

```text
python -B scripts/prepare.py import --raw-root <原始下载目录> --out <项目内新目录>
python -B scripts/prepare.py reference --measurements <measurements.csv> --tree <准确tree_id> --out <项目内新目录>
```

`--tree` 可重复，必须明确选树。只支持当前六类 IntChron JSON 和 Brehm 2021 PANGAEA 原样解压表；其他已下载格式不能假装已标准化。来源库目录结构见项目原始下载清单。

参考逐年逆方差均值；若同年离散度超过测量误差预期，标准误按 sqrt(max(1,Q/(n-1))) 扩大。该保守约定会影响结果，不能用来替代批次协方差或自动修正地区偏移。工具不会自动剔除异常值。
