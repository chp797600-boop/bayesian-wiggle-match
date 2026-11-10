# Bayesian Tree-ring Wiggle-match — 初选工具

本地MCP＋skill＋有独立日历年代依据的单年树轮¹⁴C参考。用于浮动树轮的候选年代初选与模型诊断，**不自动验收可靠绝对年代**。MCP封装0.1.0；计算接口沿用0.4.0，数学内核0.1.0。

## 内容

- `server/`：本地stdio服务、冻结计算引擎、参考数据及示例。
- `skills/bayesian-wiggle-match/`：调用指导、五关判读、偏移和失败处理。
- `server/data/catalog.json`：七份参考，合计494个年度实测点，逐树保存且附来源、年代依据、覆盖、误差限制和散列。

参考：ssm48（168点）、vahu1（107）、ua1601（51）、hib26（61）、gp22（87）、Salisbury（12）和Mickleton（8）。这不是完整年代/地区覆盖，不能把库外年份排除。默认不混合不同树或插值制造测量。示例另列，不计入494点。

## 本地运行

需要Python3.12及依赖。Windows请将仓库放在非C盘工作目录；服务输入和输出均限制在配置的工作目录内。

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r server/requirements.txt
python -B server/launch.py
```

最后一条启动stdio服务，等待MCP消息，不是交互命令行程序。客户端配置使用上述虚拟环境的Python作为`command`，`args`为`["-B", "/absolute/path/to/repository/server/launch.py"]`，`cwd`为仓库路径。默认工作目录为仓库；需要变更时通过`MCP_WORKSPACE_ROOT`指定包含服务目录的工作目录。仓库不包含作者机器的客户端配置或依赖二进制。

Codex可在可信项目的`.codex/config.toml`配置`[mcp_servers.bayesian_wiggle_match]`的command、args、cwd；按本机路径填写。skill目录复制到项目`.agents/skills/bayesian-wiggle-match`。不需要复制历史项目中的绝对路径。

## 调用流程

1. `list_reference_data`、`get_reference_info`：核查参考适用性。
2. `validate_floating_series`：检查目标CSV。
3. `run_wiggle_match`：明确参考、位置窗、偏移均值/SD及新run_id。
4. `assess_match`：按需做留一与后验预测，报告候选带和未解决项。

目标CSV列：`sample_id,ring_offset,c14_age_bp,sigma_bp,sample_span_years,measurement_ids,tree_id`。BP与正的1σ，每行一个单年轮，最早实测轮offset=0，缺测保持实际年距。天文整数年0=1 BCE。演示见`server/examples/`；浮动目标可放`inputs/`（已忽略，避免上传样品）。

偏移先验必须显式给定，示例50 BP不是通用先验。参考不足默认拒绝；显式探索时结果仅条件化于可评估位置。未知偏移的最终Acomb尚未实现；PPC不替代Acomb或独立真值检验。默认独立年度误差不处理批次协方差/年度相关，日历不确定度也未建模。

## 验证和限制

原本地项目58项测试及官方SDK完整五工具调用通过；本仓库另提供便携性检查。旧示例已知起点400 CE而MAP391，保留该不成功结果，不能用作准确性保证。不同参考、实验室、区域和生长季适用性仍需审查。仅本地Windows/stdio实际验证；其他平台需自行测试。

`skills/`内部分历史验证说明提及原项目档案；这些历史运行和论文原件不随仓库发布。参考数据的SHA和公开文献在目录中，不能将未附历史资料当作已在此仓库重跑。

## 来源与权利

参考出处：Brehm等2025（https://doi.org/10.1038/s41467-024-55757-y）、Brehm等2022（https://doi.org/10.1038/s41467-022-28804-9）、Akcicek等2026（https://doi.org/10.1111/arcm.70088）、Bayliss等2026（https://doi.org/10.1016/j.nimb.2025.165919）。逐份来源及处理记录见catalog。数据沿用各自来源许可（包括Bayliss论文CC BY-NC）；本仓库没有将第三方数据统一重新许可为MIT。源代码尚未指定开源许可。没有包含论文PDF或SDK依赖源码。
