## =============================================================================
## 更正记录（2026-10-27）：本文件原登记的 L1 发现队列为 GSE337700（骨折骨不愈计划），
## 因 NCBI 不可达未采用；实际分析使用 CELLxGene Bone Marrow Reference atlas
## （b6037d0f，正常 vs 骨关节炎）。下方 DATASETS$discovery_atlas 为实际数据源。
## =============================================================================

## =============================================================================
## <U+65B9><U+6848><U+4E09><U+FF1A><U+9053><U+5730><U+9AA8><U+4F24><U+79D1><U+836F><U+6750> + <U+9AA8><U+6298><U+4E0D><U+6108><U+5408> <U+5355><U+7EC6><U+80DE><U+865A><U+62DF><U+6572><U+9664><U+7814><U+7A76>
## <U+9879><U+76EE><U+914D><U+7F6E><U+6587><U+4EF6>
##
## <U+7528><U+9014><U+FF1A><U+96C6><U+4E2D><U+7BA1><U+7406><U+6240><U+6709><U+53EF><U+8C03><U+53C2><U+6570><U+FF0C><U+907F><U+514D><U+6563><U+843D><U+5728><U+5404><U+811A><U+672C><U+4E2D>
## <U+4FEE><U+6539><U+540E><U+9700><U+91CD><U+65B0> source <U+672C><U+6587><U+4EF6><U+624D><U+751F><U+6548>
## =============================================================================

## ---- <U+9879><U+76EE><U+8DEF><U+5F84> ----
PROJ <- list(
  root= "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa",
  data_raw        = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/data/raw",
  data_processed  = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/data/processed",
  results         = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/results",
  figures         = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/results/figures",
  tables          = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/results/tables",
  objects         = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/results/objects",
  logs            = "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/logs"
)

## ---- <U+968F><U+673A><U+79CD><U+5B50><U+FF08><U+4FDD><U+8BC1><U+53EF><U+590D><U+73B0><U+FF09> ----
SEED <- 20261008

## ---- <U+5168><U+5C40><U+5F00><U+5173> ----
## P2 <U+53EF><U+884C><U+6027><U+9A8C><U+8BC1><U+6A21><U+5F0F><U+FF1A>TRUE <U+65F6><U+53EA><U+8DD1><U+6838><U+5FC3><U+6D41><U+7A0B><U+FF0C><U+4E0D><U+505A><U+5168><U+91CF><U+5206><U+6790>
## <U+5EFA><U+8BAE><U+6B63><U+5F0F><U+5206><U+6790><U+524D><U+5148><U+7528> FALSE <U+9A8C><U+8BC1><U+6570><U+636E><U+53EF><U+884C><U+6027>
QUICK_MODE <- FALSE

## <U+865A><U+62DF><U+6572><U+9664> Tier <U+7B56><U+7565>
## Tier1 <U+7C97><U+7B5B><U+57FA><U+56E0><U+6570><U+4E0A><U+9650> / Tier2 <U+7CBE><U+7B5B><U+4E0A><U+9650> / Tier3 <U+6700><U+7EC8><U+9776><U+70B9><U+6570>
KO_TIER1_MAX <- 80L
KO_TIER2_MAX <- 20L
KO_TIER3_N   <- 5L

## ---- <U+6570><U+636E><U+96C6><U+6CE8><U+518C><U+8868> ----
## <U+6240><U+6709><U+6570><U+636E><U+96C6><U+96C6><U+4E2D><U+767B><U+8BB0><U+FF0C><U+6362><U+6570><U+636E><U+96C6><U+53EA><U+6539><U+8FD9><U+91CC>
DATASETS <- list(

  ## ============ L1 <U+53D1><U+73B0><U+5C42><U+FF1A><U+4EBA><U+9AA8><U+4E0D><U+6108> vs <U+6B63><U+5E38><U+6108><U+5408> ============
  ## GSE337700 <U+5317><U+4EAC><U+671D><U+9633><U+533B><U+9662><U+FF0C>2026-07 <U+4E0A><U+7EBF>
  ## ============ L1 发现层：人骨髓基质 atlas（实际分析使用的数据） ============
  ## 2026-10-27 更正记录：原计划用 GSE337700（骨折骨不愈，北京朝阳医院），
  ## 因 NCBI/GEO 不可达，改走 CELLxGene Discovery 的 Bone Marrow Reference
  ## atlas（b6037d0f）。实际疾病轴 = 骨关节炎(OA) vs 正常；atlas 内无骨折样本。
  discovery_atlas = list(
    portal      = "CELLxGene Discovery",
    dataset_id  = "b6037d0f-c415-45f7-8d01-0e69407fba9c",
    file        = "data/raw/bone_marrow_ref.h5ad",
    url         = "https://cellxgene.cziscience.com/",
    accessed    = "2026-10-08",
    n_cells_total = 117308L,
    disease_axis  = c("normal", "osteoarthritis"),
    normal_gse    = c("GSE190965", "GSE169396", "GSE202813"),
    oa_gse        = c("GSE147390", "GSE196678"),
    excluded_gse  = c("GSE156640", "GSE224152"),  # 血液恶性肿瘤系列，已剔除
    note = "S01-S02 由 scripts/build_s01_s02.R 复现（原始构建为交互式，见 README）"
  ),

  ## <U+4EBA><U+9AA8><U+6298> 0<U+5929> vs 30<U+5929><U+FF08><U+72EC><U+7ACB><U+4EBA> scRNA<U+FF0C><U+8DE8><U+7269><U+79CD><U+9A8C><U+8BC1><U+7528><U+FF09>
  human_timecourse_sc = list(
    gse = "GSE142786",
    note = "0<U+5929> vs 30<U+5929><U+9AA8><U+6298><U+7EC4><U+7EC7><U+FF0C><U+7528><U+4E8E><U+4EBA><U+6E90><U+8F68><U+8FF9><U+5F62><U+6001><U+9A8C><U+8BC1>"
  ),

  ## ============ L2 <U+673A><U+5236><U+5C42><U+FF1A><U+5C0F><U+9F20><U+9AA8><U+6298><U+6108><U+5408><U+65F6><U+5E8F> ============
  ## <U+65F6><U+5E8F><U+6570><U+636E><U+96C6><U+FF0C><U+63D0><U+4F9B><U+5B8C><U+6574><U+6210><U+9AA8><U+5206><U+5316><U+8F68><U+8FF9>
  mouse_timecourse_sc = list(
    ## <U+4E3B><U+65F6><U+5E8F><U+6570><U+636E><U+96C6><U+FF1A><U+4EBA><U+9AA8><U+6298> 0<U+5929> vs 30<U+5929><U+FF08>GSE142786<U+FF0C><U+5DF2><U+6838><U+5B9E><U+FF09>
    primary = list(
      gse = "GSE142786",
      species = "Homo sapiens",
      times = c(0, 30),
      verified = TRUE,
      note = "<U+5DF2><U+6838><U+5B9E><U+4E3A><U+771F><U+5B9E><U+6570><U+636E><U+96C6><U+3002><U+63D0><U+4F9B><U+5B8C><U+6574><U+6210><U+9AA8><U+5206><U+5316><U+8F68><U+8FF9>"
    ),
    ## <U+4E3B><U+65F6><U+5E8F><U+6570><U+636E><U+96C6> 2<U+FF1A><U+4EBA><U+9AA8><U+4E0D><U+6108><U+FF08>GSE337700<U+FF0C>2026-07 <U+4E0A><U+7EBF><U+FF0C><U+5317><U+4EAC><U+671D><U+9633><U+533B><U+9662><U+FF09>
    discovery = list(
      gse = "GSE337700",
      species = "Homo sapiens",
      verified = TRUE,
      note = "<U+5DF2><U+6838><U+5B9E><U+FF1A>Public on Jul 10, 2026<U+FF0C>GPL24676<U+FF0C>8 samples<U+FF0C>PRJNA1489793<U+FF0C>Citation missing"
    ),
    ## <U+5907><U+7528><U+FF1A>Broad <U+7A7A><U+95F4><U+8F6C><U+5F55><U+7EC4><U+FF08>Day0/5/15<U+FF09>
    backup_spatial = list(
      study = "SCP2948",
      times = c(0, 5, 15),
      platform = "10x Visium CytAssist",
      url = "https://singlecell.broadinstitute.org/single_cell/study/SCP2948"
    ),
    ## <U+5907><U+7528><U+FF1A>eLife Activin A <U+6570><U+636E><U+FF08>Col2-Td <U+62A5><U+544A><U+9F20><U+FF0C>d0/5/10/15/30<U+FF09>
    backup_activin = list(
      source = "eLife 89822",
      times = c(0, 5, 10, 15, 30),
      total_cells = 25429
    )
  ),

  ## ============ L3 <U+9A8C><U+8BC1><U+5C42><U+FF1A>Bulk + GWAS ============
  ## <U+9AA8><U+4E0D><U+6108> vs <U+6B63><U+5E38><U+6108><U+5408><U+FF0C><U+540C><U+5E73><U+53F0><U+6CE8><U+91CA><U+82AF><U+7247><U+FF0C><U+53EF><U+76F4><U+63A5><U+5408><U+5E76> <U+2192> <U+672C><U+9879><U+76EE><U+6700><U+91CD><U+8981><U+7684><U+72EC><U+7ACB><U+9A8C><U+8BC1><U+961F><U+5217>
  validation_bulk = list(
    nonunion = "GSE93213",
    healed   = "GSE93215",
    platform = "GPL570",           # Affymetrix HG-U133 Plus 2.0<U+FF0C><U+540C><U+5E73><U+53F0><U+53EF><U+5408><U+5E76>
    note = "Chen <U+7B49><U+539F><U+7814><U+7A76><U+4F7F><U+7528><U+3002><U+6CE8><U+610F> ANXA3 <U+5DF2><U+88AB><U+8BE5><U+7814><U+7A76><U+62A5><U+9053><U+FF0C><U+52A1><U+5FC5><U+907F><U+5F00><U+8BE5><U+57FA><U+56E0><U+4F5C><U+4E3A><U+4E3B><U+9776><U+70B9>"
  ),

  ## 2<U+578B><U+7CD6><U+5C3F><U+75C5><U+5408><U+5E76><U+9AA8><U+6298><U+FF08><U+4EE3><U+8C22><U+6027><U+9AA8><U+75C5><U+6269><U+5C55><U+961F><U+5217><U+FF09>
  dm_fracture_bulk = list(
    gse = "GSE95849",
    note = "<U+53EF><U+4F5C><U+4E3A><U+9AA8><U+4EE3><U+8C22><U+5F02><U+5E38><U+4E9A><U+578B><U+6269><U+5C55><U+5206><U+6790>"
  ),

  ## GWAS <U+6027><U+72B6><U+FF08><U+901A><U+8FC7> GWAS Catalog API <U+83B7><U+53D6><U+FF09>
  gwas_traits = c(
    "Bone mineral density",
    "Healing of fractures",
    "Bone union",
    "Osteoporosis"
  )
)

## ---- <U+649E><U+9898><U+89C4><U+907F><U+6E05><U+5355><U+FF08>2026-10-08 <U+68C0><U+7D22><U+6838><U+5B9E><U+540E><U+65B0><U+589E><U+FF0C><U+52A1><U+5FC5><U+52FF><U+6539><U+FF09> ----
## <U+4E0B><U+5217><U+57FA><U+56E0>/<U+89D2><U+5EA6><U+5DF2><U+88AB><U+5DF2><U+53D1><U+8868><U+8BBA><U+6587><U+62A2><U+5360><U+FF0C><U+51FA><U+73B0><U+5728><U+672C><U+9879><U+76EE><U+4E2D><U+4F1A><U+88AB><U+76F4><U+63A5><U+5224><U+4E3A><U+91CD><U+590D>
AVOID <- list(
  ## ---- <U+5DF2><U+88AB><U+62A2><U+5360><U+7684>"<U+79D1><U+5B66><U+95EE><U+9898><U+8868><U+8FF0>"<U+FF08><U+5199><U+4F5C><U+65F6><U+5FC5><U+987B><U+56DE><U+907F><U+8FD9><U+4E9B><U+63AA><U+8F9E><U+FF09> ----
  ## <U+2460> PMID 42353610 Curr Issues Mol Biol 2026;48(6):606<U+FF08><U+9646><U+519B><U+519B><U+533B><U+5927><U+5B66>+<U+91CD><U+5E86><U+533B><U+79D1><U+5927><U+FF09>
  ##    <U+5DF2><U+62A5><U+9053> "SSPCs arrested in an undifferentiated state in nonunion"
  ##    <U+2192> <U+672C><U+9879><U+76EE><U+7981><U+6B62><U+4F7F><U+7528>"<U+6210><U+9AA8><U+8C31><U+7CFB><U+505C><U+6EDE>/<U+672A><U+5206><U+5316>"<U+4F5C><U+4E3A><U+6838><U+5FC3><U+7ED3><U+8BBA><U+8868><U+8FF0>
  ##    <U+2192> <U+6539><U+4E3A>"<U+54EA><U+4E9B><U+57FA><U+56E0><U+7684><U+6572><U+9664><U+53EF><U+6062><U+590D><U+5206><U+5316><U+6F5C><U+80FD>"<U+FF08><U+5E72><U+9884><U+89C6><U+89D2><U+FF0C><U+975E><U+63CF><U+8FF0><U+89C6><U+89D2><U+FF09>
  ## <U+2461> PMID 40393465 Phytother Res 2025<U+FF08><U+5C71><U+4E1C><U+4E2D><U+533B><U+836F><U+5927><U+5B66><U+FF09>
  ##    <U+5DF2><U+62A5><U+9053> <U+7EED><U+65AD><U+7682><U+82F7>VI <U+2192> Piezo1 <U+2192> LEPR+ BMSC / PODXL+ EC <U+8026><U+8054>
  ##    <U+2192> <U+7981><U+6B62><U+628A> Lepr <U+4F5C><U+4E3A><U+8F68><U+8FF9><U+8D77><U+70B9><U+6807><U+8BB0><U+FF1B><U+7981><U+6B62><U+628A> Piezo1 <U+4F5C><U+4E3A><U+4E3B><U+9776><U+70B9>
  ## <U+2462> PMID 39366157 J Formul Med 2024
  ##    <U+5DF2><U+62A5><U+9053> <U+675C><U+4EF2><U+4E38> <U+2192> YAP/<U+03B2>-catenin/VEGF <U+8F74><U+FF0C>BMSC-<U+5185><U+76AE><U+5F02><U+8D28><U+6027>
  ##    <U+2192> YAP/TAZ <U+53EA><U+80FD><U+4F5C><U+4E3A><U+4E0B><U+6E38><U+80CC><U+666F><U+FF0C><U+4E0D><U+80FD><U+4F5C><U+4E3A><U+672C><U+7814><U+7A76><U+7684><U+4E3B><U+8981><U+53D1><U+73B0>
  ## <U+2463> PMID 38990653 JCI Insight 2024;9(16):e176802
  ##    <U+5DF2><U+62A5><U+9053> NF1 <U+5047><U+5173><U+8282><U+7A7A><U+95F4><U+8F6C><U+5F55><U+7EC4> + BMP <U+4FE1><U+53F7><U+53D7><U+635F>
  ##    <U+2192> <U+4E0D><U+5F97><U+518D><U+4EE5>"BMP <U+4FE1><U+53F7><U+53D7><U+635F>"<U+4F5C><U+4E3A><U+9AA8><U+4E0D><U+6108><U+7684><U+6838><U+5FC3><U+673A><U+5236><U+7ED3><U+8BBA>

  ## ---- <U+5DF2><U+88AB><U+62A2><U+5360><U+7684><U+5177><U+4F53><U+57FA><U+56E0><U+FF08><U+4E0D><U+5F97><U+4F5C><U+4E3A><U+4E3B><U+9776><U+70B9>/<U+4E3B><U+7ED3><U+8BBA><U+FF09> ----
  ## Lepr     : PMID 40393465 <U+5DF2><U+4F5C><U+4E3A><U+6838><U+5FC3><U+FF08><U+4E0E><U+8BE5><U+6587><U+91CD><U+53E0><U+5EA6><U+6700><U+9AD8><U+7684><U+8D77><U+70B9><U+57FA><U+56E0><U+FF09>
  ## Piezo1   : PMID 40393465 <U+5DF2><U+9501><U+5B9A><U+4E3A><U+836F><U+7269><U+76F4><U+63A5><U+9776><U+70B9>
  ## LIF/LIFR/IL6ST : GSE337700 <U+539F><U+6587><U+FF08><U+6570><U+636E><U+96C6><U+4F5C><U+8005><U+FF09><U+5DF2><U+62A5><U+9053><U+8BE5><U+901A><U+8BAF><U+589E><U+5F3A>
  ## Acsl4    : PMID 42353610 <U+5DF2><U+4F5C><U+4E3A><U+94C1><U+6B7B><U+4EA1><U+6838><U+5FC3>
  ## EPAS1    : PMID 42353610 <U+5DF2><U+4F5C><U+4E3A><U+4E0A><U+6E38><U+5E72><U+9884><U+9776><U+70B9>
  ## ANXA3    : GSE93213/93215 <U+539F><U+7814><U+7A76><U+5DF2><U+62A5><U+9053><U+4E3A> T2DM <U+5408><U+5E76><U+9AA8><U+6298> biomarker
  ## 人类符号（2026-10-08 修正：原为小鼠命名，在人类数据里匹不中，撞题检测会失效）
  genes = c("LEPR","PIEZO1","LIF","LIFR","IL6ST","ACSL4","EPAS1","ANXA3"),

  ## ---- <U+5199><U+4F5C><U+7981><U+7528><U+63AA><U+8F9E><U+FF08><U+6295><U+7A3F><U+524D><U+81EA><U+67E5><U+FF09> ----
  forbidden_phrases = c(
    "arrested in an undifferentiated state",
    "<U+505C><U+6EDE><U+4E8E><U+672A><U+5206><U+5316><U+72B6><U+6001>",
    "<U+8C31><U+7CFB><U+505C><U+6EDE>",
    "LEPR+ BMSC differentiation",
    "<U+9A8C><U+8BC1><U+4E86> Piezo1 <U+901A><U+8DEF>",
    "<U+9996><U+6B21>",
    "<U+56FD><U+5185><U+5916><U+9996><U+6B21>",
    "<U+586B><U+8865><U+4E86><U+7A7A><U+767D>",
    "<U+5177><U+6709><U+91CD><U+8981><U+7406><U+8BBA><U+610F><U+4E49><U+548C><U+5E94><U+7528><U+4EF7><U+503C>"
  )
)

## ---- <U+8D28><U+63A7><U+9608><U+503C> ----
QC <- list(
  ## <U+5355><U+7EC6><U+80DE>
  min_features = 200L,
  max_features = 6000L,
  min_counts   = 1000L,
  ## <U+7EBF><U+7C92><U+4F53><U+6BD4><U+4F8B><U+4E0A><U+9650><U+FF1A><U+4EBA><U+9AA8><U+7EC4><U+7EC7><U+504F><U+8106><U+FF0C><U+53D6><U+5BBD><U+677E><U+6807><U+51C6>
  max_pct_mt   = 20L,
  ## <U+6837><U+672C><U+7EA7><U+6700><U+5C0F><U+7EC6><U+80DE><U+6570><U+FF1B>GSE337700 <U+4F8B><U+5916><U+FF08><U+672C><U+5C31><U+5C11><U+FF0C><U+5FC5><U+987B><U+62A5><U+51FA><U+964D><U+7EA7><U+5904><U+7406><U+FF09>
  min_cells_per_sample = 1500L,
  small_dataset_exempt = c("NUBF02", "NUBF04", "NUBF05",
                           "G0531",  "G0601",  "G0607"),
  ## <U+6700><U+5C0F><U+8868><U+8FBE><U+6BD4><U+4F8B><U+FF08><U+57FA><U+56E0><U+5728><U+591A><U+5C11><U+7EC6><U+80DE><U+4E2D><U+8868><U+8FBE><U+624D><U+7B97><U+68C0><U+6D4B><U+5230><U+FF09>
  min_detect_pct = 0.10
)

## ---- <U+7EC6><U+80DE><U+7C7B><U+578B><U+6CE8><U+91CA><U+6807><U+8BB0><U+57FA><U+56E0> ----
## <U+6765><U+6E90><U+FF1A><U+6587><U+732E><U+5171><U+8BC6> + CellTypist/Azimuth <U+53C2><U+8003>
## =============================================================================
## 🔴🔴🔴 严重警告（2026-10-08 实测确认）🔴🔴🔴
## 下面这张 MARKERS 是【小鼠命名】（Runx2 / Col1a1 / Cd3e），
## 而本项目主数据集 bone_marrow_ref.h5ad 是【人类数据】（RUNX2 / COL1A1 / CD3E）。
##
## 实测结果（见 scripts/95_audit_markers.R）：
##   人类符号命中 9/9，小鼠符号命中 0/9
##   本表 136 个基因在人类数据里命中率 = 0.0%，24 个细胞类型全部落空！
##
## 【危害】intersect() 落空后不报错，注释静默失效 →
##   细胞全被打成 Unassigned，而流程照常跑完。这类 bug 不会崩，但结论全错。
##
## 【正确做法】人类数据用 R/02_stromal_subset.R 的 HUMAN_MARKERS。
##   用任何标记表之前，先跑 audit_markers() 自检命中率。
##
## 【保留原因】原始方案针对 GSE337700（同属人类），但此表源于小鼠文献；
##   保留并改名为 MARKERS_MOUSE + 加警告，防止误用。
## =============================================================================
MARKERS_MOUSE <- list(
  Osteoblast      = c("Runx2","Sp7","Col1a1","Bglap","Ibsp","Alpl","Phex","Dmp1","Ifitm5"),
  Osteocyte       = c("Dmp1","Phex","Ifitm5","Sost","Cox7a1a","Alpl"),
  Osteoclast      = c("Ncr1","Ctsk","Acp5","Sost","Tnfrsf11a","Dctnpp1","Atp6v0d2"),
  Chondrocyte     = c("Sox9","Acan","Col2a1","Col11a1","Hapln1","Matn3","Crim1"),
  Periosteal_Prog  = c("Prrx1","Cxcl12","Ntf3","Kcnip2","Bmp4","Col15a1"),
  Mesenchymal     = c("Vim","Pdgfra","Col1a1","Fap","Dcn","Lum"),
  Fibroblast      = c("Col1a1","Col1a2","Dcn","Lum","Pdgfra","Thbs2"),
  Endothelial     = c("Pecam1","Cdh5","Egfl7","Kdr","Vwf","Emcn"),
  Pericyte        = c("Rgs5","Pdgfrb","Cspg4","Acta2","Des"),
  Macrophage       = c("Adgre1","Cd68","Lyz2","C1qa","C1qb","Marco"),
  Monocyte        = c("Ly6c2","Ccr2","S100a8","S100a9","Plac8","Chil3"),
  T_cell          = c("Cd3e","Cd3d","Cd2","Trbc1","Lck","Thy1"),
  NK_cell         = c("Ncr1","Klrb1c","Nkg7","Prf1","Gzmb"),
  B_cell          = c("Cd79a","Cd79b","Ms4a1","Ighm","Ebf1","Pax5"),
  Plasma_cell     = c("Jchain","Mzb1","Derl3","Xbp1","Sdc1","Igkc"),
  Neutrophil       = c("S100a8","S100a9","Mpo","Elane","Camp","Lcn2"),
  HSPC            = c("Cd34","Kit","Mpl1","Gata2","Hlf","Avi1"),
  Schwann_cell    = c("Sox10","Plp1","S100b","Mpz","Pmp22"),
  Sensory_neuron  = c("Trpv1","Calca","Cgrp","Ntrk2","Ret","Scn7a"),
  NogoR_neuron    = c("Fgf9","Ngfr","Ntrk3","Runx1"),
  Smooth_muscle   = c("Acta2","Tagln","Myh11","Cnn1"),
  Tendon_cell     = c("Thbs4","Tnmd","Xirp1","Col12a1"),
  Sinusoidal_EC   = c("Lyve1","Olfact1d2","Stab2","Clec14a"),
  Arterial_EC     = c("Efnb2","Sox17","Gja5","Cxcl12")
)

## =============================================================================
## 【人类标记表】当前项目主数据集用这张（bone_marrow_ref.h5ad 为人类）
## 已实测：与数据完全匹配（人类符号命中 9/9；MARKERS_MOUSE 命中 0/9）。
## R/02_stromal_subset.R 里的 HUMAN_MARKERS 更详细（含骨髓基质专属亚群），
## 此处为全谱系通用版，供其他脚本引用。
## =============================================================================
MARKERS_HUMAN <- list(
  MSC_generic     = c("ENG","THY1","NT5E","PDGFRA","CD44","ITGB1","ALCAM"),
  Stromal_fibro   = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","FAP","THY1"),
  ## CXCL12-abundant reticular (CAR) cells —— 骨髓微环境核心支持细胞
  ## 注意：不用 LEPR 作主标记（撞题规避，见 AVOID），改用 CXCL12/ADIPOQ 轴
  CAR_cell        = c("CXCL12","ADIPOQ","KITLG","PDGFRA","VCAM1"),
  Adipo_prog      = c("ADIPOQ","PLIN1","LEP","PPARG","CEBPA","FABP4"),
  Osteo_prog      = c("RUNX2","SP7","ALPL","COL1A1","SPP1","IBSP"),
  Osteoblast      = c("RUNX2","SP7","ALPL","COL1A1","BGLAP","IBSP","SPP1","SATB2"),
  Osteocyte       = c("DMP1","PHEX","SOST","MEPE","FGF23","ENPP1"),
  Chondro_prog    = c("SOX9","ACAN","COL2A1","COL11A1","MATN3","CRTAC1"),
  Bipoten_prog    = c("PPARG","CEBPA","RUNX2","SP7","WNT5A","BMP2"),
  Endothelial     = c("PECAM1","VWF","CDH5","KDR","EGFL7","EMCN"),
  Pericyte_SMC    = c("ACTA2","TAGLN","MYH11","RGS5","PDGFRB","CNN1"),
  Hematopoietic   = c("PTPRC","CD34","KIT","MPO","ELANE"),
  Monocyte_Mac    = c("CD68","LYZ","CSF1R","ITGAM","C1QA","C1QB","CD14"),
  T_cell          = c("CD3E","CD3D","CD2","TRBC2","LCK"),
  B_cell          = c("MS4A1","CD79A","CD79B","IGHM","PAX5"),
  Plasma_cell     = c("JCHAIN","MZB1","DERL3","XBP1","SDC1"),
  NK_cell         = c("NKG7","PRF1","GZMB","KLRD1"),
  Neutrophil      = c("S100A8","S100A9","MPO","ELANE","CAMP","LCN2"),
  Erythroid       = c("HBB","HBA1","HBA2","ALAS2","GYPA")
)

## ---- 兼容别名：脚本里凡引用 MARKERS 的，一律指向人类版 ----
## 旧的 R/*.R 无需改动即可正确工作（原本引用的 MARKERS 是 0% 命中的小鼠表）。
MARKERS <- MARKERS_HUMAN

## =============================================================================
## 【污染标记表 2026-10-08 新增】
## 【坑位17】CELLxGENE 多数据集整合后 cell_type 标签不可全信：
##   实测 bone_marrow_ref.h5ad 中 cell_type=="stromal cell of bone marrow"
##   取 6379 细胞聚类后，有 1853 个（29.1%）其实是免疫/红细胞/血小板：
##     C3 浆细胞(IGHGP/MZB1)  C4 B细胞(TNFRSF13C/POU2AF1)
##     C8 髓系(ITGAX/CD86)    C10 巨噬(C1QA/B/C/MRC1)
##     C12 血小板(PF4/GP9)    C13 中性粒(MPO/ELANE/AZU1)
##     C14 T细胞(CD3E/CD3D)   C15 神经样(INSM1/ELAVL4)
##   根因：各数据集作者标注口径不一，部分把"非实质细胞"粗标为 stromal。
##   对策：不能只信 cell_type 标签，必须用 canonical marker 二次硬过滤。
## =============================================================================
CONTAMINATION <- list(
  Plasma_cell = c("JCHAIN","MZB1","DERL3","XBP1","SDC1",
                  "IGHG1","IGHA1","IGHA2","IGKC","IGLC2","IGHGP"),
  B_cell      = c("MS4A1","CD79A","CD79B","IGHM","PAX5",
                  "TNFRSF13C","POU2AF1","IRF4","SLAMF1"),
  T_cell      = c("CD3E","CD3D","CD3G","CD2","TRBC2","LCK","CD8A","IL7R","CD5"),
  Myeloid     = c("LYZ","CD68","CSF1R","ITGAM","C1QA","C1QB","C1QC",
                  "FCN1","S100A8","S100A9","AIF1","MRC1","CD14","ITGAX"),
  Neutrophil  = c("MPO","ELANE","AZU1","CAMP","PRTN3","CTSG","LTF","DEFA1"),
  Platelet    = c("PF4","PPBP","GP9","CABP5","NRGN","TMEM40"),
  Erythroid   = c("HBB","HBA1","HBA2","ALAS2","GYPA","HBD","HBM"),
  ## 泛免疫：PTPRC=CD45，所有有核血细胞的通用标记，基质细胞必须阴性
  panImmune   = c("PTPRC","CD52","LAPTM5","SRGN","CORO1A","CD37")
)

## =============================================================================
## 【分级阈值 2026-10-08 实测校准】—— 由 scripts/94_diag_contam.R 诊断得出
##
## 【坑位18】不能对所有污染类别用同一个阈值！
##   一刀切 0.25 的实测后果：33636 个基质细胞被杀掉 24070（71.6%），
##   其中 Erythroid 一项独占 13115 —— 全都是误杀。
##
## 校准方法：对每个候选阈值档，看「命中细胞」的**基质身份标记均值**。
##   基质标记仍高 → 是带背景的基质细胞（误杀）
##   基质标记≈0   → 真污染
##
## ┌──────────────┬───────────┬────────────┬────────┬──────────┐
## │ 类别         │ 0.25档命中│ 基质标记均值│ 相对全体│ 结论     │
## ├──────────────┼───────────┼────────────┼────────┼──────────┤
## │ Erythroid    │ 15992     │ 1.155      │ 1.12×   │ 完全误杀 │
## │ panImmune    │ 12715     │ 1.013      │ 0.98×   │ 误杀     │
## │ Neutrophil   │  2185     │ 0.705      │ 0.68×   │ 误杀     │
## │ Platelet     │  1745     │ 0.583      │ 0.57×   │ 可疑     │
## │ Myeloid      │  2008     │ 0.582      │ 0.56×   │ 可疑     │
## │ Plasma_cell  │  3719     │ 0.482      │ 0.47×   │ 可疑     │
## │ T_cell       │    94     │ 0.424      │ 0.41×   │ 可疑     │
## │ B_cell       │   823     │ 0.078      │ 0.08×   │ 真污染   │
## └──────────────┴───────────┴────────────┴────────┴──────────┘
##
## 【为什么血红蛋白是背景而非污染】
##   交叉表：HB 在 0.25~1 档的细胞 nFeature 中位数 1636、基质标记 1.17，
##   而 HB=0 档只有 1328 / 0.91 —— **血红蛋白高的反而是质量更好的基质细胞**。
##   骨髓样本本身含大量红细胞，裂解后游离血红蛋白 mRNA 被基质细胞吸附，
##   且细胞越大/RNA 越多吸附越多。故 Erythroid 类在本数据集中**禁用**。
##
## 【同理】中性粒标记 S100A8/S100A9 在基质细胞中有真实基础表达，须放宽。
## =============================================================================
CONTAMINATION_THRESHOLD <- list(
  Plasma_cell = 0.5,    # 0.5 档：2320 命中，基质标记 0.199（0.19×）→ 真污染
  B_cell      = 0.25,   # 823 命中，基质标记 0.078（0.08×）→ 很确定
  T_cell      = 0.5,    # 命中极少（94），取 0.5 保守
  Myeloid     = 0.5,    # 0.5 档：863 命中，基质标记 0.260（0.25×）→ 真污染
  Neutrophil  = 1.0,    # 各档基质标记都 >0.5，放宽到 1.0（仅 104 命中）
  Platelet    = 1.0,    # 1.0 档：150 命中，基质标记 0.227（0.22×）→ 真污染
  Erythroid   = Inf,    # ★ 禁用：骨髓基质普遍背景，见上方说明
  panImmune   = 1.0     # 1.0 档：533 命中，基质标记 0.235（0.23×）→ 真污染
)

## ---- 基质身份标记：真基质细胞必表达，用于「双条件」判定 ----
## 真污染细胞必须同时满足：谱系标记高 **且** 基质身份标记低
STROMA_IDENTITY <- c("COL1A1","COL1A2","COL3A1","DCN","LUM","CXCL12",
                     "PDGFRA","THY1","ENG","NT5E","ACTA2","MYH11",
                     "SP7","RUNX2","ALPL")

## ---- <U+62DF><U+65F6><U+5E8F><U+8F68><U+8FF9><U+914D><U+7F6E> ----
## <U+649E><U+9898><U+4FEE><U+6B63><U+FF08>2026-10-08<U+FF09><U+FF1A><U+5DF2><U+79FB><U+9664> Lepr
## <U+7406><U+7531><U+FF1A>PMID 40393465<U+FF08>Phytother Res 2025<U+FF09><U+5DF2><U+628A>"LEPR+ BMSC <U+8F68><U+8FF9>"<U+4F5C><U+4E3A><U+6838><U+5FC3><U+6846><U+67B6><U+FF0C>
##       <U+82E5><U+672C><U+9879><U+76EE><U+518D><U+4EE5> LEPR <U+4E3A><U+8D77><U+70B9><U+FF0C><U+7B49><U+4E8E><U+4E0E><U+8BE5><U+6587><U+6B63><U+9762><U+540C><U+6784>
## <U+6539><U+7528> Cxcl12 <U+4F5C><U+4E3A><U+9759><U+6B62><U+671F><U+9AA8><U+9AD3><U+57FA><U+8D28><U+7EC6><U+80DE><U+6807><U+8BB0><U+FF08>PMID 31949165 <U+5DF2><U+8BC1> CXCL12+ BMSC
##     <U+5728><U+635F><U+4F24><U+540E><U+624D><U+8F6C><U+5206><U+5316><U+FF0C><U+53EF><U+4F5C><U+4E3A><U+8D77><U+70B9><U+FF0C><U+4F46><U+672C><U+6587><U+4E0D><U+505A><U+8F68><U+8FF9><U+4E3B><U+7EBF><U+FF0C><U+4EC5><U+4F5C><U+7EC6><U+80DE><U+5B9A><U+4E49><U+FF09>
## <U+4FEE><U+6B63><U+540E><U+8D77><U+70B9><U+6807><U+8BB0><U+FF1A>Prrx1<U+FF08><U+9AA8><U+819C><U+6765><U+6E90><U+FF09>+ Cxcl12<U+FF08><U+9759><U+6B62> BMSC<U+FF09>+ Ntf3/Kcnip2<U+FF08><U+4FDD><U+7559><U+FF09>
TRAJECTORY <- list(
  ## <U+8D77><U+70B9><U+FF1A><U+9AA8><U+819C><U+6765><U+6E90><U+6210><U+9AA8><U+524D><U+4F53> + <U+9759><U+6B62><U+671F><U+9AA8><U+9AD3><U+57FA><U+8D28><U+7EC6><U+80DE>
  ## <U+6CE8><U+610F><U+FF1A>Lepr <U+5DF2><U+88AB><U+79FB><U+51FA><U+FF08><U+649E><U+9898><U+FF09><U+FF0C><U+8BE6><U+89C1> AVOID <U+6CE8><U+91CA>
  ## 人类符号（2026-10-08 修正：原为小鼠命名 Prrx1/Cxcl12/Ntf3/Kcnip2，
  ## 在人类数据中命中率 0%，会导致轨迹起点判定完全失效）
  ## 【2026-10-08 二次修正：按实测 Spearman 相关性重选】
  ##   实测（n=6155，成骨轨迹）：
  ##     CXCL12 rho=-0.355 阳性率 92.4%  → 随拟时序下降 ✓
  ##     LPL    rho=-0.228 阳性率 61.3%  → 下降 ✓
  ##     KITLG  rho=-0.184 阳性率 56.3%  → 下降 ✓
  ##   但原来写的 PRRX1/PDGFRA/THY1 实测是**上升**的
  ##     （PRRX1 +0.304、PDGFRA +0.221、THY1 +0.237）——
  ##     因为轨迹中段是 Stromal_fibro（成骨祖细胞），这三个基因在中段达峰，
  ##     不是起点标记。NTF3 阳性率 0.0%、KCNIP2 仅 1.2%，信号太弱不可用。
  root_markers = c("CXCL12","KITLG","LPL"),
  ## <U+5DF2><U+5254><U+9664><U+7684><U+649E><U+9898><U+57FA><U+56E0>
  removed_root_markers = c("LEPR"),

  ## 【end_markers 为什么必须换掉 DMP1/PHEX/SOST】
  ##   原写 DMP1/PHEX/IFITM5/SOST/ALPL，实测全部不合格：
  ##     DMP1 rho=+0.075 阳性率 0.2%   ← 骨细胞标记
  ##     PHEX rho=-0.061 阳性率 1.0%   ← 骨细胞标记
  ##     SOST rho=+0.030 阳性率 0.0%   ← 骨细胞标记
  ##     ALPL rho=-0.186 阳性率 50.8%  ← CAR 细胞里就高表达，非成骨特异
  ##   **骨髓基质里没有骨细胞**（骨细胞包埋在骨基质中，骨髓抽吸取不到），
  ##   拿骨细胞标记做"成骨轨迹终点"是方法论错误，必然相关不上。
  ##
  ## 换成实测验证过的成骨细胞/成骨基质标记（8/8 全部强正相关）：
  ##     WIF1 +0.437  IBSP +0.424  OMD +0.405  SFRP4 +0.381
  ##     BGLAP +0.366  SPP1 +0.341  CHAD +0.296  IFITM5 +0.294
  end_markers  = c("IBSP","BGLAP","SPP1","IFITM5","OMD","WIF1","SFRP4","CHAD"),
  ## 中段成骨祖细胞标记（PDGFRA/THY1 在中段达峰，见上方说明）
  mid_markers  = c("RUNX2","SP7","COL1A1","PDGFRA","THY1"),
  ## DDRTree <U+53C2><U+6570><U+FF08><U+7981><U+7528> UMAP/MDS <U+505A><U+8F68><U+8FF9><U+FF0C><U+4F1A><U+626D><U+66F2><U+62D3><U+6251><U+FF09>
  ntops = 200L,
  norm_gamma = 2,
  use_partition = FALSE,
  ## <U+5206><U+5C94><U+70B9><U+68C0><U+6D4B>
  branch_genes_min_entropy = 20L,
  ## Monocle3 <U+5206><U+8FA8><U+7387><U+53C2><U+6570>
  learn_graph_control = 100,
  use_partition_root_leaf = FALSE
)

## ---- <U+865A><U+62DF><U+6572><U+9664><U+914D><U+7F6E> ----
KO_CONFIG <- list(
  ## scTenifoldKnk <U+53C2><U+6570>
  n_cells_per_sample = 500L,   # scTenifoldNet <U+62BD><U+6837><U+7EC6><U+80DE><U+6570>
  n_neighbors = 30L,
  nc_n_comp = 3L,
  td_K = 3L,
  td_maxIter = 1000L,
  td_nDecimal = 3L,
  td_max_error = 0.1,
  td_max_iter = 1000L,
  ## <U+5DEE><U+5F02><U+8C03><U+63A7><U+57FA><U+56E0><U+7B5B><U+9009>
  ko_fdr = 0.05,
  ## <U+53EF><U+5E72><U+9884><U+6027><U+8BC4><U+5206><U+6743><U+91CD>
  steer_weight_gene  = 0.4,    # KO <U+6548><U+5E94><U+57FA><U+56E0><U+4E0E><U+5F02><U+5E38><U+57FA><U+56E0><U+91CD><U+53E0><U+5EA6>
  steer_weight_fate  = 0.4,    # CellRank <U+547D><U+8FD0><U+6F5C><U+80FD><U+6062><U+590D><U+5EA6>
  steer_weight_track = 0.2,    # <U+8F68><U+8FF9><U+4F4D><U+7F6E><U+56DE><U+9000><U+5EA6>
  ## <U+5019><U+9009><U+9776><U+70B9><U+7EB3><U+5165><U+6807><U+51C6>
  min_steer_score = 0.5,
  ## Augur <U+654F><U+611F><U+6027><U+FF08><U+72EC><U+7ACB><U+7B2C><U+4E09><U+65B9><U+5DE5><U+5177><U+FF09>
  run_augur = TRUE,
  augur_top_n = 50L
)

## ---- <U+7F51><U+7EDC><U+836F><U+7406><U+5B66><U+914D><U+7F6E> ----
NP_CONFIG <- list(
  ## AD <U+7B5B><U+9009><U+9608><U+503C><U+FF08>TCMSP <U+7ECF><U+5178><U+9608><U+503C><U+FF09>
  ad_absorb_threshold  = 4.5,   # <U+5438><U+6536><U+5EA6> >= 4.5 (LogP)
  ad_drug_likeness     = 0.3,   # <U+836F><U+7269>likeness >= 0.3
  ad_ob_threshold      = 30,    # <U+53E3><U+670D><U+751F><U+7269><U+5229><U+7528><U+5EA6> >= 30%
  ad_caco2_threshold   = -5.15, # Caco-2 <U+6E17><U+900F><U+6027> LogP >= -5.15
  ad_bbb_threshold     = -0.3,   # BBB <U+6E17><U+900F><U+6027> >= -0.3
  ## <U+6570><U+636E><U+5E93><U+7248><U+672C><U+FF08><U+8BBA><U+6587><U+4E2D><U+5FC5><U+987B><U+5199><U+660E><U+FF09>
  db_versions = list(
    TCMSP = "v2.1",
    ETCM  = "v2.0",
    HERB  = "2023",
    SwissTargetPrediction = "2024-02"
  )
)

## ---- <U+76EE><U+6807><U+836F><U+6750><U+4E0E><U+6210><U+5206><U+FF08><U+5BD2><U+5730><U+9053><U+5730><U+9AA8><U+4F24><U+79D1><U+836F><U+6750><U+FF09> ----
## <U+6CE8><U+610F><U+FF1A><U+5DF2><U+5254><U+9664><U+56DE><U+9633><U+7389><U+9F99><U+818F><U+FF08><U+542B><U+9E9D><U+9999><U+FF0C><U+6D89><U+53CA><U+4FDD><U+62A4><U+52A8><U+7269><U+4F26><U+7406><U+4E89><U+8BAE><U+FF09>
##
## <U+649E><U+9898><U+4FEE><U+6B63><U+FF08>2026-10-08 <U+68C0><U+7D22><U+6838><U+5B9E><U+FF09><U+FF1A>
##   <U+9AA8><U+788E><U+8865><U+FF08><U+7EED><U+65AD><U+FF09><U+4ECE>"<U+5E76><U+5217><U+4E3B><U+4F53>"<U+964D><U+7EA7><U+4E3A>"<U+5BF9><U+7167><U+836F><U+6750>"
##   <U+7406><U+7531><U+FF1A>PMID 40393465<U+FF08>Phytother Res 2025<U+FF09><U+5DF2><U+7528><U+7EED><U+65AD><U+7682><U+82F7>VI <U+505A><U+8FC7>
##         "LEPR+ BMSC + <U+5206><U+5B50><U+5BF9><U+63A5> + Piezo1"<U+FF0C><U+4E0E><U+9AA8><U+788E><U+8865><U+8DEF><U+7EBF><U+9AD8><U+5EA6><U+91CD><U+5408>
##   <U+523A><U+4E94><U+52A0><U+5347><U+4E3A><U+552F><U+4E00><U+4E3B><U+6253><U+836F><U+6750>
##   <U+7406><U+7531><U+FF1A><U+68C0><U+7D22><U+672A><U+89C1>"<U+523A><U+4E94><U+52A0> + <U+9AA8><U+6298><U+6108><U+5408> + <U+5355><U+7EC6><U+80DE>"<U+7684><U+4EFB><U+4F55><U+5148><U+4F8B>
HERBS <- list(
  ciwujia = list(          # <U+523A><U+4E94><U+52A0> <U+2014> <U+5BD2><U+5730><U+9F99><U+836F><U+4E4B><U+738B><U+FF0C><U+672C><U+9879><U+76EE><U+552F><U+4E00><U+4E3B><U+6253>
    cn = "<U+523A><U+4E94><U+52A0>",
    latin = "Eleutherococcus senticosus",
    role = "primary",       # <U+4E3B><U+6253>
    key_compounds = c(
      "<U+523A><U+4E94><U+52A0><U+82F7>E", "<U+7D2B><U+4E01><U+9999><U+82F7>", "<U+5F02><U+55EA><U+76AE><U+5576>",
      "<U+7EFF><U+539F><U+9178>", "<U+91D1><U+7EBF><U+83B2><U+82F7>", "<U+69F2><U+76AE><U+7D20>"
    ),
    novelty = "<U+68C0><U+7D22><U+672A><U+53D1><U+73B0><U+523A><U+4E94><U+52A0> + <U+9AA8><U+6298><U+6108><U+5408> scRNA-seq <U+7684><U+5148><U+4F8B>"
  ),
  gushubi = list(          # <U+9AA8><U+788E><U+8865> <U+2014> <U+964D><U+7EA7><U+4E3A><U+5BF9><U+7167><U+836F><U+6750>
    cn = "<U+9AA8><U+788E><U+8865>",
    latin = "Drynaria roosii",
    role = "comparator",    # <U+5BF9><U+7167><U+FF08><U+649E><U+9898><U+964D><U+7EA7><U+FF09>
    key_compounds = c("<U+85AF><U+84E3><U+7682><U+82F7>", "<U+69F2><U+76AE><U+7D20>", "<U+5C71><U+67F0><U+915A>", "<U+7EFF><U+539F><U+9178>", "<U+963F><U+9B4F><U+9178>"),
    collision_note = "PMID 40393465 <U+5DF2><U+7528><U+7EED><U+65AD><U+7682><U+82F7>VI <U+8D70>LEPR+/Piezo1 <U+8DEF><U+7EBF><U+FF0C><U+672C><U+9879><U+76EE><U+4EC5><U+4F5C><U+5BF9><U+7167>"
  ),
  lurong = list(           # <U+9E7F><U+8338> <U+2014> <U+586B><U+7CBE><U+8865><U+9AD3><U+FF0C><U+80BE><U+4E3B><U+9AA8><U+751F><U+9AD3>
    cn = "<U+9E7F><U+8338>",
    latin = "P Cervi Cornu Pantophysi",
    key_compounds = c("<U+9E7F><U+8338><U+591A><U+80BD>", "<U+9E7F><U+89D2><U+80F6>", "<U+603B><U+86CB><U+767D>", "<U+7F9F><U+812F><U+6C28><U+9178>"),
    note = "<U+591A><U+80BD><U+7C7B><U+6210><U+5206><U+9700><U+5355><U+72EC><U+5904><U+7406><U+5206><U+5B50><U+5BF9><U+63A5><U+FF08><U+4E0D><U+9002><U+5408><U+6807><U+51C6><U+5C0F><U+5206><U+5B50> docking<U+FF09>"
  ),
  zuoguiwan = list(        # <U+5DE6><U+5F52><U+4E38> <U+2014> <U+80BE><U+7CBE><U+4E8F><U+865A><U+8BC1><U+FF08><U+9AA8><U+4E0D><U+6108><U+5408><U+6838><U+5FC3><U+8BC1><U+578B><U+FF09>
    cn = "<U+5DE6><U+5F52><U+4E38>",
    latin = "Zuo Gui Wan",
    type = "formula",
    note = "<U+65B9><U+5242><U+7EA7><U+FF0C><U+9700><U+4ECE> TCMSP <U+65B9><U+5242><U+5E93><U+83B7><U+53D6><U+6210><U+5206><U+96C6><U+5408>"
  ),
  buyanghuanwu = list(     # <U+8865><U+9633><U+8FD8><U+4E94><U+6C64> <U+2014> <U+6C14><U+865A><U+8840><U+7600><U+578B>
    cn = "<U+8865><U+9633><U+8FD8><U+4E94><U+6C64>",
    latin = "Bu Yang Huan Wu Tang",
    type = "formula",
    note = "<U+65B9><U+5242><U+7EA7><U+FF0C><U+8865><U+6C14><U+6D3B><U+8840><U+901A><U+7EDC>"
  )
)

## ---- <U+5168><U+5C40><U+4E3B><U+9898><U+FF08>ggplot<U+FF09> ----
THEME <- list(
  base_size = 11,
  colors = list(
    nonunion = "#A32D2D",       # <U+9AA8><U+4E0D><U+6108> - <U+7EA2>
    healed   = "#1D9E75",       # <U+6B63><U+5E38><U+6108><U+5408> - <U+7EFF>
    osteoblast = "#185FA5",
    chondrocyte = "#BA7517",
    osteoclast  = "#993C1D",
    fibroblast  = "#5F5E5A",
    macrophage  = "#D4537E",
    endothelial = "#0F6E56",
    neuron      = "#534AB7",
    highlight   = "#BA7517"
  )
)

## ---- <U+5206><U+6790><U+9636><U+6BB5><U+5F00><U+5173> ----
STEPS <- list(
  s01_qc          = TRUE,
  s02_integrate   = TRUE,
  s03_annotate    = TRUE,
  s04_trajectory  = TRUE,
  s05_regnetwork  = TRUE,
  s06_ko_screen    = TRUE,
  s07_ko_fate     = TRUE,
  s08_network_pharm = TRUE,
  s09_docking     = TRUE,
  s10_validate    = TRUE,
  s11_figures     = TRUE
)

message("<U+914D><U+7F6E><U+6587><U+4EF6><U+5DF2><U+52A0><U+8F7D><U+FF1A>",
        PROJ$root,
        " | SEED=", SEED,
        " | QUICK_MODE=", QUICK_MODE)
