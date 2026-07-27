# te-core — 生产性能核心（Rust）

de novo 重复检测的性能关键路径,按 V4_DESIGN §3 用 Rust 实现(原型阶段的 Python 只做编排,本核心替换其中不 scale 的热路径)。分阶段构建,每阶段可独立编译 + 对照验证。

## Phase 1 ✅ `te-seed` — io + canonical k-mer 索引 + 种子采集

替换 `../src/proto_scaled.py` 的 Stage 0–2(numpy 扫基因组 + jellyfish + occurrence 收集)。

- **io**:FASTA 流式读、按记录边界隔离(k-mer 不跨记录)、64-bit 坐标。
- **index**:canonical (k≤16) k-mer,增量 fwd/rc 滚动编码;每 k-mer 打包为 `(canon:u32<<32)|gpos:u32` 单 u64,std unstable sort 分组 → 计数与采集为一次线性扫描(无外部 crate,离线可编译;rayon 并行留 Phase 3)。
- **seed**:count≥t 的 repetitive 种子 + **hash-bottom-k 出现次数封顶**(确定性、顺序无关)。

**验证(TAIR10, k=16):** 六项统计与 jellyfish **完全一致**——Total 133,916,481 / Distinct 92,998,560 / Unique 79,662,342 / count=2 9,163,114 / R(≥3) 4,173,104 / Max 47,238;count≥200 种子 6,983(同 jellyfish);occurrence 位置抽验 3/3 正确(含 ± 链);封顶 ≤200/种子。**运行时 5.6 s、峰值 1.18 GB**(io+index+seed 一趟完成)。

```
cargo build --release
./target/release/te-seed <genome.fa> [k=16] [min_count=3] [cap=200] [--seeds out.tsv]
```

## Phase 2 ✅ `te-discover` — A1 extend-to-growing-consensus(替换 cd-hit)

cd-hit all-vs-all 聚类是多 Gb 的悬崖(本项目实测撞过:~99 万 instances grinding;minimap2 all-vs-all 在重复核上爆炸)。A1 用 RepeatScout/RECON 范式 O(occ) 取代:高频种子(贪心,count 降序)→ 每 occurrence 只对**左右生长的 consensus 比对一次**(逐列多数 fit),持续匹配的=家族成员、发散的丢弃;家族的基因组跨度被 mask,后续种子不重建。**这正是 cd-hit 该换的那块;POA consensus 仍交 spoa(见下 Phase 2b)。**

- **成员招募 = 聚类**:extend 中持续 mismatch 的 occurrence 被踢出(非本家族),存活的就是家族——无两两比较。
- **§E6 串联守卫**:`is_tandem`(周期 2–600 自相似 >60%)把串联/卫星分流到 `*.tandem.fa`,不进 interspersed 库(否则 12 kb runaway consensi 会拖垮 RepeatMasker——本项目实测验证)。

**验证(TAIR10, count≥200):** 667 家族 → **182 interspersed + 485 tandem(分流)**(T2T 完整着丝粒 → 高频种子多为卫星,符合预期)。**231,869 次成员→consensus 比对(O(occ),非两两)**,7 s / 1.2 GB。
**mask 率:mdl-repeat + A1 = 24.17%**,对照 cd-hit+spoa 原型 pass1 的 24.32%——**差 ~0.15 pp,A1 内核作为 cd-hit 替代验证通过**。

```
./target/release/te-discover <genome.fa> [k=16] [min_count=200] [cap=200] [out_prefix]
# -> <prefix>.consensi.fa (interspersed) + <prefix>.tandem.fa + <prefix>.members.bed
```

> 诚实边界:A1 consensus 是 **ungapped 粗 consensus**(182 中 12 个 blastn 全长)——mask 率达标是因 RepeatMasker 用粗 consensus 也能贴回拷贝;但全长一致序列质量需 **Phase 2b** 把 `members.bed` 的成员跨度交 **spoa** 做 gapped POA 精修(即「spoa 留」那一半,预期把全长率从 12 拉高、如原型 Stage4 的 5/40→21/33)。

## Phase 2b ✅ `te-refine` — spoa gapped POA 精修

「在成熟工具失败处自建(A1 聚类),在它们擅长处复用(POA)」:A1 输出 `members.bed`,本步逐家族提取成员序列(按链定向)交 **spoa** 做 gapped POA,得高质量一致序列(ungapped A1 consensus 比不了的)。Rust 编排,spoa 做活,不重写 POA。

**验证(TAIR10):** A1 ungapped(a1c)24.17% → **A1 + spoa 精修(a1r)24.26%**(+0.09 pp),对照 cd-hit+spoa 原型 24.32%——**差 0.06 pp**。即 **Rust A1(替 cd-hit)+ spoa(POA)复现原型 mask 率**,整套「自建+复用」架构验证通过。

```
./target/release/te-refine <genome.fa> <members.bed> <out.fa>   # 调 spoa -r 0 逐家族
```
> 注:`is_tandem` 只查 consensus **内部**周期,**大周期卫星**(consensus=单个卫星单元、无内部重复)会漏判进 interspersed 库——守卫待补(加「单元在基因组内串联拷贝数」判据)。不影响 mask 率,只影响 interspersed 纯度。

## Phase 3 ✅ rayon 并行 + minimizer 规模化

- **rayon 并行**:按记录并行扫描 + `par_sort_unstable`。**exact(w=1)仍完全对齐 jellyfish**(Distinct 92,998,560 / R 4,173,104),**3.95 s vs 单线程 5.6 s**。(并行 concat 峰值内存升到 ~2.7 GB,多 Gb 靠 minimizer 解。)
- **minimizer(`--w`,单调 deque,1/w 密度)= 多 Gb 内存/吞吐杠杆**:w=11 → 24.3 M minimizers(1/6 密度)、**1.33 s / 589 MB**(vs exact 2.7 GB);高拷贝家族种子作为 minimizer 被保留。discovery 仍工作(w=11:66 interspersed + 377 tandem,70k 比对,vs exact 182+485/232k——更少家族、~3× 更省,w 可调密度/召回)。

```
./target/release/te-seed     <genome.fa> 16 3 200 --w 11   # minimizer 模式
./target/release/te-discover <genome.fa> 16 200 200 out --w 11
```

**多 Gb 缩放实测(1.34 Gb 合成基因组 = 10× TAIR10,转 Rust 的核心理由):**

| | time | peak RAM | → 15 Gb 外推 |
|---|---|---|---|
| exact (w=1) | 56 s | 40.5 GB | ~450 GB(悬崖) |
| **minimizer (w=11)** | 11 s | **8.4 GB** | **~94 GB(可行)** |

exact-mode 内存 ∝ 基因组长度 → 15 Gb 撞墙;**minimizer 把它压成可行(~5× 省内存、~5× 更快)**,验证了多 Gb 可扩展性这一转 Rust 的初衷。

**te-discover(A1 引擎,非仅索引)在 1.34 Gb 实测:** w=11 下 **~90 s / 5.4 GB**——贪心 A1 循环 + masked 向量在多 Gb 上内存有界(→ 15 Gb 约 ~60 GB,可行)。

**高拷贝家族「mask 全拷贝」修复:** A1 原先只 mask 建库采样的 200 个成员 → 高拷贝家族(如 Alu ~10⁶ 拷贝)被反复重建。修复:建 consensus 后用 `member_match`(对 consensus ≥70% ungapped 比对)**招募并 mask 全部真拷贝**(不碰无关 seed 命中,避免过度 mask)。TAIR10 家族数从「采样-mask 157(冗余)/全-seed-mask 15(过度)」收敛到 **92(recruit-mask,正解)**,mask 率 23.93%(vs 采样-mask 24.17%、原型 24.32%)——少 40% 家族、mask 仅 −0.24 pp,冗余大降而召回基本保持。

**大周期卫星守卫(已补):** `is_tandem` 只查 consensus 内部周期,漏判单元卫星;新增 `genomic_tandem`——成员跨度在基因组里 >50% 串联相邻 → 卫星,分流 tandem 轨。实测把 25 个单元卫星从 interspersed 移出(182→157),分散的高拷贝元件正确保留。

## `dtr` — Pan_TE 集成 CLI(取代 repgraph/dtr 占位)

Pan_TE 早有 te-looker 槽位(`dtr run --genome <fa> --out <dir> --threads <n>` → `<dir>/families.fasta`),但 `dtr` 二进制从未建出。`dtr` 二进制把已验证的组件**串联**起来填这个槽:

```
dtr run --genome <fa> --out <dir> [--threads N] [--min-count C] [--w W]
        [--dfam-hmm H --enable-known-hmm-discovery] [--final-gated]
  └─ te-discover (A1 extend-to-consensus,替 cd-hit) → disc.consensi + disc.members.bed
  └─ te-refine   (spoa gapped POA 精修 de novo members) → families.fasta
  └─ dtr-hmm     (Track 1: nhmmer --tblout → overlap/copy gate → known_hmm.members.bed)
  └─ te-refine   (Track 1 POA) → known_hmm.families.fasta → append 至 raw families
  └─ --final-gated: dtr-stitch → dtr-gate/CD-HIT → families.fasta + te-looker-gated-v1 provenance

Track 2 partial tools: dtr-track2 wires dtr-window → dtr-lsh → dtr-community → dtr-component-poa → dtr-nonte-guard(optional, manifest-auditable) → dtr-copy → dtr-structure → dtr-accept-audit; all outputs retain family_call=false / mode=partial, acceptance_complete=false, and final_library=null.
```

Track 1 只在同时给出 `--dfam-hmm` 与 `--enable-known-hmm-discovery` 时运行。`dtr-hmm` 以 0-based half-open BED 输出命中；先在同一模型/染色体的重叠域中保留最高分命中（不将双链同位点重复计为两个拷贝），再要求 ≥5 个长度 ≥80 bp 的成员。最终仍由统一 copy/length gate 与跨 seed merge 决定是否提升到输出库。

**接进 Pan_TE(已做):** ① `Pan_TE/bin/` 软链 `dtr`/`te-discover`/`te-refine`/`te-seed` → `core/target/release/*`;② `Pan_TE/bin/Pan_TE` 的 `resolve_te_looker_bin` 加了 `BIN_DIR/dtr` 兜底,使其像其它检测器一样自动从 `bin/` 找到 dtr。下游 `combine_results` 的 `cd-hit-est` 正常消费 `families.fasta`。
> 前置:`cd core && cargo build --release`; Track 1 需要 `nhmmer`，精修需要 `spoa`。HMM Track-1 已实现；graph window-stride/Track 2 仍为后续工作。

## 后续(规划)

- Track 2 已有 H0 spaced-seed 成本门禁、全局 hash-bottom 候选窗口、受限 MinHash/LSH、确定性连通容器、Louvain-style 局部社区细分、组件 POA、BLASTN copy catalog 及显式 non-TE reference guard。拟南芥 rRNA guard 在 7 条 provisional consensus 中拒绝了 6 条，说明 copy/POA 不足以证明 TE。所有中间物均明确 `mode=partial` / `family_call=false`，绝不进入最终库；非-TE panel 已可覆盖细胞器、rRNA、tRNA 与 snRNA；蛋白/宿主基因守卫仍是完整接纳的缺口。下一阶段是完整 Leiden aggregation、由目前的保守结构审计驱动的边界重估、LTR/TIR/Helitron 的专用证据及 G1–G6 审计，之后才可接入最终 family gate。其后为把 A1/refine 重构进 lib 让 `dtr` 全程 in-process(免 subprocess)、sweep-line 覆盖记账、striped-lock 哈希、exact concat 零拷贝、真实多 Gb 端到端 mask 率。

设计依据:`../docs/V4_DESIGN.md` §A/§3;原型与 masking 实证:`../docs/DESIGN_REVIEW.md`。
