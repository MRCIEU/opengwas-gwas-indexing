SET GLOBAL max_connections = 500;

ALTER INSTANCE DISABLE INNODB REDO_LOG;

create table if not exists opengwas.phewas
(
    id        bigint unsigned auto_increment,
    gwas_id_n mediumint unsigned null,
    snp_id    varchar(255)       null,
    chr_id    tinyint unsigned   not null,
    pos       int unsigned       not null,
    ea        varchar(255)       null,
    nea       varchar(255)       null,
    eaf       float              null,
    beta      float              null,
    se        float              null,
    lp        float              null comment '-log10(p), coerced to 999999 when p was provided as inf i.e. when -log10(p) was provided as string 0',
    ss        varchar(255)       null,
    primary key (id, chr_id, pos)
)
    partition by range (((`chr_id` * 1000000000) + `pos`)) (

      partition p_chr1_01 values less than (1014000000),
      partition p_chr1_02 values less than (1027000000),
      partition p_chr1_03 values less than (1042000000),
      partition p_chr1_04 values less than (1056000000),
      partition p_chr1_05 values less than (1070000000),
      partition p_chr1_06 values less than (1084000000),
      partition p_chr1_07 values less than (1099000000),
      partition p_chr1_08 values less than (1112000000),
      partition p_chr1_09 values less than (1154000000),
      partition p_chr1_10 values less than (1167000000),
      partition p_chr1_11 values less than (1181000000),
      partition p_chr1_12 values less than (1195000000),
      partition p_chr1_13 values less than (1209000000),
      partition p_chr1_14 values less than (1222000000),
      partition p_chr1_15 values less than (1236000000),
      partition p_chr1_16 values less than (2000000000),

      partition p_chr2_01 values less than (2013000000),
      partition p_chr2_02 values less than (2026000000),
      partition p_chr2_03 values less than (2038000000),
      partition p_chr2_04 values less than (2050000000),
      partition p_chr2_05 values less than (2062000000),
      partition p_chr2_06 values less than (2076000000),
      partition p_chr2_07 values less than (2098000000),
      partition p_chr2_08 values less than (2113000000),
      partition p_chr2_09 values less than (2128000000),
      partition p_chr2_10 values less than (2143000000),
      partition p_chr2_11 values less than (2159000000),
      partition p_chr2_12 values less than (2173000000),
      partition p_chr2_13 values less than (2187000000),
      partition p_chr2_14 values less than (2203000000),
      partition p_chr2_15 values less than (2217000000),
      partition p_chr2_16 values less than (2231000000),
      partition p_chr2_17 values less than (3000000000),

      partition p_chr3_01 values less than (3012000000),
      partition p_chr3_02 values less than (3025000000),
      partition p_chr3_03 values less than (3039000000),
      partition p_chr3_04 values less than (3053000000),
      partition p_chr3_05 values less than (3066000000),
      partition p_chr3_06 values less than (3081000000),
      partition p_chr3_07 values less than (3098000000),
      partition p_chr3_08 values less than (3112000000),
      partition p_chr3_09 values less than (3126000000),
      partition p_chr3_10 values less than (3141000000),
      partition p_chr3_11 values less than (3155000000),
      partition p_chr3_12 values less than (3169000000),
      partition p_chr3_13 values less than (3185000000),
      partition p_chr3_14 values less than (4000000000),

      partition p_chr4_01 values less than (4011000000),
      partition p_chr4_02 values less than (4024000000),
      partition p_chr4_03 values less than (4037000000),
      partition p_chr4_04 values less than (4054000000),
      partition p_chr4_05 values less than (4067000000),
      partition p_chr4_06 values less than (4080000000),
      partition p_chr4_07 values less than (4094000000),
      partition p_chr4_08 values less than (4107000000),
      partition p_chr4_09 values less than (4121000000),
      partition p_chr4_10 values less than (4135000000),
      partition p_chr4_11 values less than (4150000000),
      partition p_chr4_12 values less than (4164000000),
      partition p_chr4_13 values less than (4179000000),
      partition p_chr4_14 values less than (5000000000),

      partition p_chr5_01 values less than (5014000000),
      partition p_chr5_02 values less than (5030000000),
      partition p_chr5_03 values less than (5044000000),
      partition p_chr5_04 values less than (5061000000),
      partition p_chr5_05 values less than (5078000000),
      partition p_chr5_06 values less than (5094000000),
      partition p_chr5_07 values less than (5109000000),
      partition p_chr5_08 values less than (5121000000),
      partition p_chr5_09 values less than (5135000000),
      partition p_chr5_10 values less than (5152000000),
      partition p_chr5_11 values less than (5166000000),
      partition p_chr5_12 values less than (6000000000),

      partition p_chr6_01 values less than (6013000000),
      partition p_chr6_02 values less than (6026000000),
      partition p_chr6_03 values less than (6030000000),
      partition p_chr6_04 values less than (6032000000),
      partition p_chr6_05 values less than (6033000000),
      partition p_chr6_06 values less than (6041000000),
      partition p_chr6_07 values less than (6055000000),
      partition p_chr6_08 values less than (6071000000),
      partition p_chr6_09 values less than (6086000000),
      partition p_chr6_10 values less than (6100000000),
      partition p_chr6_11 values less than (6115000000),
      partition p_chr6_12 values less than (6130000000),
      partition p_chr6_13 values less than (6145000000),
      partition p_chr6_14 values less than (6159000000),
      partition p_chr6_15 values less than (7000000000),

      partition p_chr7_01 values less than (7010000000),
      partition p_chr7_02 values less than (7022000000),
      partition p_chr7_03 values less than (7036000000),
      partition p_chr7_04 values less than (7050000000),
      partition p_chr7_05 values less than (7067000000),
      partition p_chr7_06 values less than (7082000000),
      partition p_chr7_07 values less than (7098000000),
      partition p_chr7_08 values less than (7113000000),
      partition p_chr7_09 values less than (7130000000),
      partition p_chr7_10 values less than (7146000000),
      partition p_chr7_11 values less than (8000000000),

      partition p_chr8_01 values less than (8008000000),
      partition p_chr8_02 values less than (8014000000),
      partition p_chr8_03 values less than (8025000000),
      partition p_chr8_04 values less than (8040000000),
      partition p_chr8_05 values less than (8058000000),
      partition p_chr8_06 values less than (8073000000),
      partition p_chr8_07 values less than (8088000000),
      partition p_chr8_08 values less than (8104000000),
      partition p_chr8_09 values less than (8119000000),
      partition p_chr8_10 values less than (8133000000),
      partition p_chr8_11 values less than (9000000000),

      partition p_chr9_01 values less than (9011000000),
      partition p_chr9_02 values less than (9024000000),
      partition p_chr9_03 values less than (9037000000),
      partition p_chr9_04 values less than (9084000000),
      partition p_chr9_05 values less than (9098000000),
      partition p_chr9_06 values less than (9113000000),
      partition p_chr9_07 values less than (9127000000),
      partition p_chr9_08 values less than (10000000000),

      partition p_chr10_01 values less than (10012000000),
      partition p_chr10_02 values less than (10024000000),
      partition p_chr10_03 values less than (10038000000),
      partition p_chr10_04 values less than (10056000000),
      partition p_chr10_05 values less than (10068000000),
      partition p_chr10_06 values less than (10082000000),
      partition p_chr10_07 values less than (10096000000),
      partition p_chr10_08 values less than (10109000000),
      partition p_chr10_09 values less than (10123000000),
      partition p_chr10_10 values less than (11000000000),

      partition p_chr11_01 values less than (11011000000),
      partition p_chr11_02 values less than (11025000000),
      partition p_chr11_03 values less than (11039000000),
      partition p_chr11_04 values less than (11051000000),
      partition p_chr11_05 values less than (11066000000),
      partition p_chr11_06 values less than (11080000000),
      partition p_chr11_07 values less than (11094000000),
      partition p_chr11_08 values less than (11108000000),
      partition p_chr11_09 values less than (11121000000),
      partition p_chr11_10 values less than (12000000000),

      partition p_chr12_01 values less than (12015000000),
      partition p_chr12_02 values less than (12029000000),
      partition p_chr12_03 values less than (12045000000),
      partition p_chr12_04 values less than (12060000000),
      partition p_chr12_05 values less than (12076000000),
      partition p_chr12_06 values less than (12093000000),
      partition p_chr12_07 values less than (12108000000),
      partition p_chr12_08 values less than (12121000000),
      partition p_chr12_09 values less than (13000000000),

      partition p_chr13_01 values less than (13032000000),
      partition p_chr13_02 values less than (13046000000),
      partition p_chr13_03 values less than (13060000000),
      partition p_chr13_04 values less than (13074000000),
      partition p_chr13_05 values less than (13088000000),
      partition p_chr13_06 values less than (13102000000),
      partition p_chr13_07 values less than (14000000000),

      partition p_chr14_01 values less than (14035000000),
      partition p_chr14_02 values less than (14050000000),
      partition p_chr14_03 values less than (14065000000),
      partition p_chr14_04 values less than (14079000000),
      partition p_chr14_05 values less than (14094000000),
      partition p_chr14_06 values less than (15000000000),

      partition p_chr15_01 values less than (15038000000),
      partition p_chr15_02 values less than (15052000000),
      partition p_chr15_03 values less than (15063000000),
      partition p_chr15_04 values less than (15077000000),
      partition p_chr15_05 values less than (15090000000),
      partition p_chr15_06 values less than (16000000000),

      partition p_chr16_01 values less than (16010000000),
      partition p_chr16_02 values less than (16025000000),
      partition p_chr16_03 values less than (16056000000),
      partition p_chr16_04 values less than (16070000000),
      partition p_chr16_05 values less than (16081000000),
      partition p_chr16_06 values less than (17000000000),

      partition p_chr17_01 values less than (17013000000),
      partition p_chr17_02 values less than (17031000000),
      partition p_chr17_03 values less than (17044000000),
      partition p_chr17_04 values less than (17053000000),
      partition p_chr17_05 values less than (17068000000),
      partition p_chr17_06 values less than (18000000000),

      partition p_chr18_01 values less than (18014000000),
      partition p_chr18_02 values less than (18035000000),
      partition p_chr18_03 values less than (18050000000),
      partition p_chr18_04 values less than (18064000000),
      partition p_chr18_05 values less than (19000000000),

      partition p_chr19_01 values less than (19011000000),
      partition p_chr19_02 values less than (19022000000),
      partition p_chr19_03 values less than (19037000000),
      partition p_chr19_04 values less than (19048000000),
      partition p_chr19_05 values less than (20000000000),

      partition p_chr20_01 values less than (20015000000),
      partition p_chr20_02 values less than (20034000000),
      partition p_chr20_03 values less than (20049000000),
      partition p_chr20_04 values less than (21000000000),

      partition p_chr21_01 values less than (21032000000),
      partition p_chr21_02 values less than (22000000000),

      partition p_chr22_01 values less than (22030000000),
      partition p_chr22_02 values less than (22041000000),
      partition p_chr22_03 values less than (23000000000),

      partition p_chrX_01 values less than (23038000000),
      partition p_chrX_02 values less than (23079000000),
      partition p_chrX_03 values less than (23117000000),
      partition p_chrX_04 values less than (24000000000),

      partition p_chrY_01 values less than (25000000000),

      partition p_chrMT_01 values less than (26000000000),

      partition p_chrMax_01 values less than (MAXVALUE)
    );

create index phewas_idx_gwas_id_n
    on opengwas.phewas (gwas_id_n);

-- create index phewas_idx_cpalleles
--     on opengwas.phewas (chr_id, pos, ea, nea);


create table opengwas.`tophits_5e-8_10000_0.001`
(
    id        bigint unsigned auto_increment,
    gwas_id_n mediumint unsigned null,
    snp_id    varchar(255)       null,
    chr_id    tinyint unsigned   not null,
    pos       int unsigned       not null,
    ea        varchar(255)       null,
    nea       varchar(255)       null,
    eaf       float              null,
    beta      float              null,
    se        float              null,
    lp        float              null comment '-log10(p), coerced to 999999 when p was provided as inf i.e. when -log10(p) was provided as string 0',
    ss        varchar(255)       null,
    primary key (id)
);

create index tophits_idx_gwas_id_n
    on opengwas.`tophits_5e-8_10000_0.001` (gwas_id_n);

create table opengwas.`tophits_1e-5_1000_0.8`
(
    id        bigint unsigned auto_increment,
    gwas_id_n mediumint unsigned null,
    snp_id    varchar(255)       null,
    chr_id    tinyint unsigned   not null,
    pos       int unsigned       not null,
    ea        varchar(255)       null,
    nea       varchar(255)       null,
    eaf       float              null,
    beta      float              null,
    se        float              null,
    lp        float              null comment '-log10(p), coerced to 999999 when p was provided as inf i.e. when -log10(p) was provided as string 0',
    ss        varchar(255)       null,
    primary key (id)
);

create index tophits_idx_gwas_id_n
    on opengwas.`tophits_1e-5_1000_0.8` (gwas_id_n);


-- Check number of records in total in phewas
SELECT
    sum(TABLE_ROWS)
FROM
    information_schema.PARTITIONS
WHERE
    TABLE_SCHEMA = 'opengwas'
    AND TABLE_NAME = 'phewas';


-- Check number of records by partitions in phewas
SELECT
    PARTITION_NAME,
    TABLE_ROWS
FROM
    information_schema.PARTITIONS
WHERE
    TABLE_SCHEMA = 'opengwas'
    AND TABLE_NAME = 'phewas';


-- Check number of records by chr in phewas
SELECT
    SUBSTRING_INDEX(PARTITION_NAME, '_', 2) AS partition_prefix,
    SUM(TABLE_ROWS) AS total_rows
FROM
    information_schema.PARTITIONS
WHERE
    TABLE_SCHEMA = 'opengwas'
    AND TABLE_NAME = 'phewas'
GROUP BY
    partition_prefix
ORDER BY
    partition_prefix;


-- Importing:
-- Upgrade service instance to 24C 96G
-- Change resources limits to 22C 90G in docker compose yml
-- Upgrade OCI Block Volume to Higher performance (VPU/GB:20) i.e. 50k IOPS, 614 MB/s
-- Takes 2 hours

-- Indexing:
-- Downgrade service instance to 8C 32G
-- Upgrade OCI Block Volume to Balanced performance (VPU/GB:10) i.e. 25k IOPS, 480 MB/s
-- Takes 6+ hours


PURGE BINARY LOGS TO '';

ALTER INSTANCE ENABLE INNODB REDO_LOG;
