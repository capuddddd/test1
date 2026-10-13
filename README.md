# Single-cell RNA analysis notes

GSE121638 单细胞 RNA 分析笔记和图像。入口为 [analysis.md](analysis.md)，
该文档包含需要按步骤执行的 R 代码和示例输出，不是一个完整自动化流水线。

## 输入整理

```bash
python3 scripts/organize_10x.py ./data/GSE121638_RAW
python3 scripts/organize_10x.py ./data/GSE121638_RAW --apply
```

默认预览；`--apply` 才移动文件。脚本跳过目录、符号链接和不匹配的文件，
按 GU 数字加 P/T 的样本名分组；目标已存在时拒绝覆盖。
执行移动前应有数据备份，输入和样本名需与 GEO 数据说明核对。

## 分析环境

笔记涉及 Seurat、harmony、clustree、scCustomize、ggsci 等 R 包；
需从真实成功运行的环境记录版本和 sessionInfo()。
QC、批次校正和 marker 注释应按实验设计核实。
图像使用仓库相对路径，支持本地查看。

本仓库未包含完整依赖锁、数据访问/使用说明和 LICENSE；
作者需要确认这些信息后再宣称可复现或完成许可审核。

## 本次修复及验证边界

修复 KEGG 部分多余右括号、S4/data.frame 混用、缺失 ratio helper 和缺少的包导入。
热图改用现有 Seurat 对象及 marker，替换未定义的 `st.data` / `markGenes`；
使用副本保留细胞类型标签，Slingshot 图标签按实际簇数循环。
PC 选择保留原启发式规则，遇到无候选值明确报错；细胞注释长度与簇数不符也会报错。
这些修复没有通过真实 R 数据运行验证，原图不代表修改后代码的新结果。

```bash
python3 -m unittest discover -s tests -v
```
