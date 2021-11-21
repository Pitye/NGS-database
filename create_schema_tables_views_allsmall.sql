CREATE SCHEMA ngs_forensic DEFAULT CHARACTER SET latin1 ;

CREATE TABLE ngs_forensic.nomenclature_autostr (
  locus varchar(100) NOT NULL,
  PubMed_ID varchar(100) NOT NULL,
  seq_name varchar(500) NOT NULL,
  allele varchar(100) NOT NULL,
  PubMed_Nomenclature_AutoSTR_Illumina_longseq varchar(500) NOT NULL,
  reverse_Illumina_longseq varchar(500) DEFAULT NULL,
  Illumina_sequence varchar(500) NOT NULL,
  present_in_Illumina varchar(100) NOT NULL,
  PRIMARY KEY (PubMed_ID)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

create table ngs_forensic.heads(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL unique,
project varchar(50) NOT NULL,
analysis varchar(50) NOT NULL,
run varchar(50) NOT NULL,
gender varchar(50) NOT NULL,
created varchar(100) NOT NULL,                             
no_mismatches int,
PRIMARY KEY (id)
);

create table ngs_forensic.heads_family_tree(
  id int NOT NULL AUTO_INCREMENT,
  sample_name varchar(50) NOT NULL,
  project varchar(50) NOT NULL,
  analysis varchar(50) NOT NULL,
  run varchar(50) NOT NULL,
  gender varchar(50) NOT NULL,
  created varchar(100) NOT NULL,
  no_mismatches int DEFAULT NULL,
  PRIMARY KEY (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

create table ngs_forensic.heads_y(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL unique,
project varchar(50) NOT NULL,
analysis varchar(50) NOT NULL,
run varchar(50) NOT NULL,
gender varchar(50) NOT NULL,
created varchar(100) NOT NULL,                             
no_mismatches_y int,
PRIMARY KEY (id)
);

create table ngs_forensic.heads_x(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL unique,
project varchar(50) NOT NULL,
analysis varchar(50) NOT NULL,
run varchar(50) NOT NULL,
gender varchar(50) NOT NULL,
created varchar(100) NOT NULL,                             
no_mismatches_x int,
PRIMARY KEY (id)
);

CREATE TABLE ngs_forensic.heads_flankingreg (
  id int NOT NULL AUTO_INCREMENT,
  sample_name varchar(50) NOT NULL,
  project varchar(50) NOT NULL,
  analysis varchar(50) NOT NULL,
  run varchar(50) NOT NULL,
  gender varchar(50) NOT NULL,
  created varchar(100) NOT NULL,
  no_mismatches int DEFAULT NULL,
  PRIMARY KEY (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE ngs_forensic.heads_y_flankingreg (
  id int NOT NULL AUTO_INCREMENT,
  sample_name varchar(50) NOT NULL,
  project varchar(50) NOT NULL,
  analysis varchar(50) NOT NULL,
  run varchar(50) NOT NULL,
  gender varchar(50) NOT NULL,
  created varchar(100) NOT NULL,
  no_mismatches_y int DEFAULT NULL,
  PRIMARY KEY (id),
  UNIQUE KEY sample_name (sample_name)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

create table ngs_forensic.autostrdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
sequence varchar(500),
CE_validation varchar(50),
head_id int,
PRIMARY KEY (id),
CONSTRAINT head_fk FOREIGN KEY (head_id) REFERENCES heads(id) ON DELETE CASCADE
);

create table ngs_forensic.y_strdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
sequence varchar(500),
CE_validation varchar(50),
head_y_id int,
PRIMARY KEY (id),
CONSTRAINT head_y_fk FOREIGN KEY (head_y_id) REFERENCES heads_y(id) ON DELETE CASCADE
);

create table ngs_forensic.x_strdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
sequence varchar(500),
CE_validation varchar(50),
head_x_id int,
PRIMARY KEY (id),
CONSTRAINT head_x_fk FOREIGN KEY (head_x_id) REFERENCES heads_x(id) ON DELETE CASCADE
);

create table ngs_forensic.snpdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
validation_by_no_reads varchar(50),
head_id int,
PRIMARY KEY (id),
CONSTRAINT head_snp_fk FOREIGN KEY (head_id) REFERENCES heads(id) ON DELETE CASCADE
);

CREATE TABLE ngs_forensic.autostrdata_flankingreg (
  id int NOT NULL AUTO_INCREMENT,
  sample_name varchar(50) NOT NULL,
  marker varchar(50) NOT NULL,
  allele varchar(50) DEFAULT NULL,
  no_reads int DEFAULT NULL,
  sequence varchar(500) DEFAULT NULL,
  CE_validation varchar(50) DEFAULT NULL,
  head_id int DEFAULT NULL,
  PRIMARY KEY (id),
  KEY head_fk_FR (head_id),
  CONSTRAINT head_fk_FR FOREIGN KEY (head_id) REFERENCES heads_flankingreg (id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE ngs_forensic.autostrdata_family_tree (
  id int NOT NULL AUTO_INCREMENT,
  sample_name varchar(50) NOT NULL,
  marker varchar(50) NOT NULL,
  allele varchar(50) DEFAULT NULL,
  no_reads int DEFAULT NULL,
  sequence varchar(500) DEFAULT NULL,
  CE_validation varchar(50) DEFAULT NULL,
  head_id int DEFAULT NULL,
  PRIMARY KEY (id),
  KEY head_fk_FT (head_id),
  CONSTRAINT head_fk_FT FOREIGN KEY (head_id) REFERENCES heads_family_tree (id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE ngs_forensic.y_strdata_flankingreg (
  id int NOT NULL AUTO_INCREMENT,
  sample_name varchar(50) NOT NULL,
  marker varchar(50) NOT NULL,
  allele varchar(50) DEFAULT NULL,
  no_reads int DEFAULT NULL,
  sequence varchar(500) DEFAULT NULL,
  CE_validation varchar(50) DEFAULT NULL,
  head_y_id int DEFAULT NULL,
  PRIMARY KEY (id),
  KEY head_y_fk_FR (head_y_id),
  CONSTRAINT head_y_fk_FR FOREIGN KEY (head_y_id) REFERENCES heads_y_flankingreg (id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

create view ngs_forensic.markerautostrview as 
select m.marker, m.allele,  m.sequence,  m. avg_no_reads, m.count_seq, m.count_seq/n.sum_count_seq  as frequency 
from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from autostrdata AS a
			JOIN ( select sequence, marker,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from autostrdata
                group by  sequence, marker) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker 
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as m
join  (select marker,  sum(count_seq) as sum_count_seq
			from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from autostrdata AS a
			JOIN ( select sequence, marker,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from autostrdata
                group by  sequence, marker) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as x
            group by marker) as n
on m.marker = n.marker 
order by marker, allele ; 

CREATE VIEW ngs_forensic.markerautostrview_flankingreg 
AS 
  SELECT m.marker, 
         m.allele, 
         m.SEQUENCE, 
         m. avg_no_reads, 
         m.count_seq, 
         m.count_seq / n.sum_count_seq AS frequency 
  FROM   (SELECT DISTINCT ( a.SEQUENCE ), 
                          a.allele, 
                          a.marker, 
                          t.count_seq, 
                          t.avg_no_reads 
          FROM   autostrdata_flankingreg AS a 
                 join (SELECT SEQUENCE, 
                              marker, 
                              Count(SEQUENCE) AS count_seq, 
                              Avg(no_reads)   AS avg_no_reads 
                       FROM   autostrdata_flankingreg 
                       GROUP  BY SEQUENCE, 
                                 marker) AS t 
                   ON a.SEQUENCE = t.SEQUENCE 
                      AND a.marker = t.marker 
          WHERE  a.ce_validation = 'validated_CE' 
                  OR a.ce_validation = 'validated_no_reads') AS m 
         join (SELECT marker, 
                      SUM(count_seq) AS sum_count_seq 
               FROM   (SELECT DISTINCT ( a.SEQUENCE ), 
                                       a.allele, 
                                       a.marker, 
                                       t.count_seq, 
                                       t.avg_no_reads 
                       FROM   autostrdata_flankingreg AS a 
                              join (SELECT SEQUENCE, 
                                           marker, 
                                           Count(SEQUENCE) AS count_seq, 
                                           Avg(no_reads)   AS avg_no_reads 
                                    FROM   autostrdata_flankingreg 
                                    GROUP  BY SEQUENCE, 
                                              marker) AS t 
                                ON a.SEQUENCE = t.SEQUENCE 
                                   AND a.marker = t.marker 
                       WHERE  a.ce_validation = 'validated_CE' 
                               OR a.ce_validation = 'validated_no_reads') AS x 
               GROUP  BY marker) AS n 
           ON m.marker = n.marker 
  ORDER  BY marker, 
            allele; 

create view ngs_forensic.markery_strview as 
select m.marker, m.allele,  m.sequence,  m. avg_no_reads, m.count_seq, m.count_seq/n.sum_count_seq  as frequency 
from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from y_strdata AS a
			JOIN ( select sequence, marker,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from y_strdata
                group by  sequence, marker) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker 
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as m
join  (select marker,  sum(count_seq) as sum_count_seq
			from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from y_strdata AS a
			JOIN ( select sequence, marker, 
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from y_strdata
                group by  sequence, marker) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker 
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as x
            group by marker) as n
on m.marker = n.marker 
order by marker, allele ;  


create view ngs_forensic.markerx_strview as 
select m.marker, m.allele,  m.sequence,  m. avg_no_reads, m.count_seq, m.count_seq/n.sum_count_seq  as frequency 
from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from x_strdata AS a
			JOIN ( select sequence, marker, 
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from x_strdata
                group by  sequence, marker) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker 
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as m
join  (select marker,  sum(count_seq) as sum_count_seq
			from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from x_strdata AS a
			JOIN ( select sequence, marker,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from x_strdata
                group by  sequence, marker) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker 
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as x
            group by marker) as n
on m.marker = n.marker 
order by marker, allele ;  

create view ngs_forensic.markersnpview as 
select t1.marker, t1.allele,  t1.avg_no_reads, t1.count_allele, t1.count_allele/t2.sum_count_allele as frequency
from (select t.allele, t.marker, t.validation_by_no_reads, count(t.allele) AS count_allele, avg(t.no_reads) AS avg_no_reads
			from snpdata as t
			where validation_by_no_reads = 'validated_no_reads'
			group by  allele, marker, validation_by_no_reads
			) as t1
join ( select t1.marker, t1.validation_by_no_reads, sum(t1.count_allele) as sum_count_allele
		from (select t.allele, t.marker, t.validation_by_no_reads, count(t.allele) AS count_allele, avg(t.no_reads) AS avg_no_reads
			from snpdata as t
			where validation_by_no_reads = 'validated_no_reads'
			group by  allele, marker, validation_by_no_reads
			) as t1
			group by  marker, validation_by_no_reads ) as t2
 on  t1.marker = t2.marker
 order by marker, allele; 

CREATE VIEW ngs_forensic.nomenclature_autostrdata AS
    SELECT 
        ASD.sample_name AS sample_name,
        ASD.marker AS marker,
        ASD.allele AS allele,
        NAY.seq_name AS seq_name,
        NAY.PubMed_ID AS PubMed_ID,
        ASD.no_reads AS no_reads,
        ASD.CE_validation AS CE_validation,
        ASD.head_id AS head_id,
        ASD.sequence AS sequence
    FROM
        (ngs_forensic.autostrdata ASD
        LEFT JOIN (SELECT 
            NA.locus AS locus,
                NA.PubMed_ID AS PubMed_ID,
                NA.seq_name AS seq_name,
                NA.PubMed_Nomenclature_AutoSTR_Illumina_longseq AS PubMed_Nomenclature_AutoSTR_Illumina_longseq,
                NA.reverse_Illumina_longseq AS reverse_Illumina_longseq,
                NA.Illumina_sequence AS Illumina_sequence
        FROM
            ngs_forensic.nomenclature_autostr NA
        WHERE
            (NA.present_in_Illumina = 'yes')) NAY ON (((ASD.sequence = NAY.Illumina_sequence)
            AND (ASD.marker = NAY.locus))));
            
            

CREATE VIEW ngs_forensic.nomenclature_markerautostrview AS
    SELECT 
        m.marker AS marker,
        m.allele AS allele,
        n.seq_name AS seq_name,
        m.avg_no_reads AS avg_no_reads,
        m.count_seq AS count_seq,
        m.frequency AS frequency,
        m.sequence AS sequence,
        n.PubMed_Nomenclature_AutoSTR_Illumina_longseq AS PubMed_Nomenclature_AutoSTR_Illumina_longseq,
        n.PubMed_ID AS PubMed_ID
    FROM
        (ngs_forensic.markerautostrview m
        LEFT JOIN (SELECT 
            n1.locus AS locus,
                n1.PubMed_ID AS PubMed_ID,
                n1.seq_name AS seq_name,
                n1.allele AS allele,
                n1.PubMed_Nomenclature_AutoSTR_Illumina_longseq AS PubMed_Nomenclature_AutoSTR_Illumina_longseq,
                n1.reverse_Illumina_longseq AS reverse_Illumina_longseq,
                n1.Illumina_sequence AS Illumina_sequence,
                n1.present_in_Illumina AS present_in_Illumina
        FROM
            ngs_forensic.nomenclature_autostr n1
        WHERE
            (n1.present_in_Illumina = 'yes')) n ON (((m.sequence = n.Illumina_sequence)
            AND (m.marker = n.locus))))
    ORDER BY m.marker , m.allele;
		 
CREATE VIEW ngs_forensic.nomenclature_markerautostrview_flankingreg 
AS 
  SELECT t.marker                                       AS marker, 
         t.allele                                       AS allele, 
         t.seq_name                                     AS seq_name, 
         t.avg_no_reads                                 AS avg_no_reads, 
         t.count_seq                                    AS count_seq, 
         t.frequency                                    AS frequency, 
         t.sequence                                     AS sequence, 
         t.pubmed_id                                    AS PubMed_ID, 
         t.pubmed_nomenclature_autostr_illumina_longseq AS 
         PubMed_Nomenclature_AutoSTR_Illumina_longseq 
  FROM   (SELECT m.marker                                       AS marker, 
                 m.allele                                       AS allele, 
                 n.seq_name                                     AS seq_name, 
                 m.avg_no_reads                                 AS avg_no_reads, 
                 m.count_seq                                    AS count_seq, 
                 m.frequency                                    AS frequency, 
                 m.sequence                                     AS sequence, 
                 n.pubmed_nomenclature_autostr_illumina_longseq AS 
                 PubMed_Nomenclature_AutoSTR_Illumina_longseq, 
                 n.pubmed_id                                    AS PubMed_ID, 
                 n.present_in_illumina                          AS 
                 present_in_Illumina 
          FROM   markerautostrview_flankingreg AS m 
                 LEFT JOIN nomenclature_autostr AS n 
                        ON n.locus = m.marker 
                           AND n.allele = m.allele 
          WHERE  n.pubmed_nomenclature_autostr_illumina_longseq = m.sequence 
                  OR n.reverse_illumina_longseq = m.sequence) AS t 
  WHERE  t.present_in_illumina = 'yes' 
          OR t.present_in_illumina = 'no'; 
  
    
CREATE VIEW ngs_forensic.nomen_freq_autostrdata AS
    SELECT 
        b.sample_name AS sample_name,
        b.marker AS marker,
        b.allele AS allele,
        b.seq_name AS seq_name,
        b.PubMed_ID AS PubMed_ID,
        b.sequence AS sequence,
        b.no_reads AS no_reads,
        b.CE_validation AS CE_validation,
        b.head_id AS head_id,
        c.avg_no_reads AS avg_no_reads,
        c.count_seq AS count_seq,
        c.frequency AS frequency
    FROM
        (ngs_forensic.nomenclature_autostrdata b
        LEFT JOIN ngs_forensic.markerautostrview c ON (((c.sequence = b.sequence)
            AND (c.marker = b.marker))))
    ORDER BY b.sample_name;

 CREATE VIEW ngs_forensic.nomenclature_autostrdata_flankingreg 
AS 
  SELECT t.sample_name   AS sample_name, 
         t.marker        AS marker, 
         t.allele        AS allele, 
         t.seq_name      AS seq_name, 
         t.pubmed_id     AS PubMed_ID, 
         t.no_reads      AS no_reads, 
         t.ce_validation AS CE_validation, 
         t.head_id       AS head_id, 
         t.sequence      AS sequence 
  FROM   (SELECT a.sample_name         AS sample_name, 
                 a.marker              AS marker, 
                 a.allele              AS allele, 
                 n.seq_name            AS seq_name, 
                 n.pubmed_id           AS PubMed_ID, 
                 a.no_reads            AS no_reads, 
                 a.ce_validation       AS CE_validation, 
                 a.head_id             AS head_id, 
                 a.sequence            AS sequence, 
                 n.present_in_illumina AS present_in_Illumina 
          FROM   ngs_forensic.nomenclature_autostr AS n 
                 LEFT JOIN autostrdata_flankingreg AS a 
                        ON n.locus = a.marker 
                           AND n.allele = a.allele 
          WHERE  n.pubmed_nomenclature_autostr_illumina_longseq = a.sequence 
                  OR n.reverse_illumina_longseq = a.sequence) AS t 
  WHERE  t.present_in_illumina = 'yes' 
          OR t.present_in_illumina = 'no';  
          
 CREATE VIEW ngs_forensic.nomenclature_autostrdata_family_tree 
AS 
  SELECT t.sample_name   AS sample_name, 
         t.marker        AS marker, 
         t.allele        AS allele, 
         t.seq_name      AS seq_name, 
         t.pubmed_id     AS PubMed_ID, 
         t.no_reads      AS no_reads, 
         t.ce_validation AS CE_validation, 
         t.head_id       AS head_id, 
         t.sequence      AS sequence 
  FROM   (SELECT a.sample_name         AS sample_name, 
                 a.marker              AS marker, 
                 a.allele              AS allele, 
                 n.seq_name            AS seq_name, 
                 n.pubmed_id           AS PubMed_ID, 
                 a.no_reads            AS no_reads, 
                 a.ce_validation       AS CE_validation, 
                 a.head_id             AS head_id, 
                 a.sequence            AS sequence, 
                 n.present_in_illumina AS present_in_Illumina 
          FROM   ngs_forensic.nomenclature_autostr AS n 
                 LEFT JOIN autostrdata_family_tree AS a 
                        ON n.locus = a.marker 
                           AND n.allele = a.allele 
          WHERE  n.pubmed_nomenclature_autostr_illumina_longseq = a.sequence 
                  OR n.reverse_illumina_longseq = a.sequence) AS t 
  WHERE  t.present_in_illumina = 'yes' 
          OR t.present_in_illumina = 'no'; 
 
 CREATE VIEW ngs_forensic.nomen_freq_autostrdata_flankingreg AS
    SELECT 
        b.sample_name AS sample_name,
        b.marker AS marker,
        b.allele AS allele,
        b.seq_name AS seq_name,
        b.PubMed_ID AS PubMed_ID,
        b.sequence AS sequence,
        b.no_reads AS no_reads,
        b.CE_validation AS CE_validation,
        b.head_id AS head_id,
        c.avg_no_reads AS avg_no_reads,
        c.count_seq AS count_seq,
        c.frequency AS frequency
    FROM
        (ngs_forensic.nomenclature_autostrdata_flankingreg b
        LEFT JOIN ngs_forensic.markerautostrview_flankingreg c ON (((c.sequence = b.sequence)
            AND (c.marker = b.marker))))
    ORDER BY b.sample_name;
 
  CREATE VIEW ngs_forensic.nomen_freq_autostrdata_family_tree AS
    SELECT 
        b.sample_name AS sample_name,
        b.marker AS marker,
        b.allele AS allele,
        b.seq_name AS seq_name,
        b.PubMed_ID AS PubMed_ID,
        b.sequence AS sequence,
        b.no_reads AS no_reads,
        b.CE_validation AS CE_validation,
        b.head_id AS head_id,
        c.avg_no_reads AS avg_no_reads,
        c.count_seq AS count_seq,
        c.frequency AS frequency
    FROM
        (ngs_forensic.nomenclature_autostrdata_family_tree b
        LEFT JOIN ngs_forensic.markerautostrview_flankingreg c ON (((c.sequence = b.sequence)
            AND (c.marker = b.marker))))
    ORDER BY b.sample_name;
    
 CREATE VIEW ngs_forensic.noname_sequence AS
    SELECT 
        nm.marker AS marker,
        nm.allele AS allele,
        ad.CE_validation AS CE_validation,
        nm.sequence AS sequence,
        nm.avg_no_reads AS avg_no_reads,
        nm.count_seq AS count_seq,
        nm.frequency AS frequency,
        ad.sample_name AS sample_name
    FROM
        (((SELECT 
            a.marker AS marker,
                a.allele AS allele,
                a.seq_name AS seq_name,
                a.avg_no_reads AS avg_no_reads,
                a.count_seq AS count_seq,
                a.frequency AS frequency,
                a.sequence AS sequence,
                a.PubMed_Nomenclature_AutoSTR_Illumina_longseq AS PubMed_Nomenclature_AutoSTR_Illumina_longseq,
                a.PubMed_ID AS PubMed_ID
        FROM
            ngs_forensic.nomenclature_markerautostrview a
        WHERE
            ISNULL(a.seq_name))) nm
        JOIN ngs_forensic.autostrdata ad ON (((nm.sequence = ad.sequence)
            AND (nm.marker = ad.marker))))
    ORDER BY nm.marker , nm.allele;
    
              

create view ngs_forensic.markery_strview_flankingreg 
AS 
  SELECT m.marker, 
         m.allele, 
         m.SEQUENCE, 
         m. avg_no_reads, 
         m.count_seq, 
         m.count_seq / n.sum_count_seq AS frequency 
  FROM   (SELECT DISTINCT ( a.SEQUENCE ), 
                          a.allele, 
                          a.marker, 
                          t.count_seq, 
                          t.avg_no_reads 
          FROM   y_strdata_flankingreg AS a 
                 join (SELECT SEQUENCE, 
                              marker, 
                              Count(SEQUENCE) AS count_seq, 
                              Avg(no_reads)   AS avg_no_reads 
                       FROM   y_strdata_flankingreg 
                       GROUP  BY SEQUENCE, 
                                 marker) AS t 
                   ON a.SEQUENCE = t.SEQUENCE 
                      AND a.marker = t.marker 
          WHERE  a.ce_validation = 'validated_CE' 
                  OR a.ce_validation = 'validated_no_reads') AS m 
         join (SELECT marker, 
                      SUM(count_seq) AS sum_count_seq 
               FROM   (SELECT DISTINCT ( a.SEQUENCE ), 
                                       a.allele, 
                                       a.marker, 
                                       t.count_seq, 
                                       t.avg_no_reads 
                       FROM   y_strdata_flankingreg AS a 
                              join (SELECT SEQUENCE, 
                                           marker, 
                                           Count(SEQUENCE) AS count_seq, 
                                           Avg(no_reads)   AS avg_no_reads 
                                    FROM   y_strdata_flankingreg 
                                    GROUP  BY SEQUENCE, 
                                              marker) AS t 
                                ON a.SEQUENCE = t.SEQUENCE 
                                   AND a.marker = t.marker 
                       WHERE  a.ce_validation = 'validated_CE' 
                               OR a.ce_validation = 'validated_no_reads') AS x 
               GROUP  BY marker) AS n 
           ON m.marker = n.marker 
  ORDER  BY marker, 
            allele; 



CREATE VIEW ngs_forensic.freq_y_strdata_flankingreg AS
    SELECT 
        b.sample_name AS sample_name,
        b.marker AS marker,
        b.allele AS allele,
        b.sequence AS sequence,
        b.no_reads AS no_reads,
        b.CE_validation AS CE_validation,
        b.head_y_id AS head_id,
        c.avg_no_reads AS avg_no_reads,
        c.count_seq AS count_seq,
        c.frequency AS frequency
    FROM
        (ngs_forensic.markery_strview_flankingreg c
        LEFT JOIN ngs_forensic.y_strdata_flankingreg b ON (((c.sequence = b.sequence)
            AND (c.marker = b.marker))))
    ORDER BY b.sample_name;

--- create table ngs_forensic.markersautostr(
id INT NOT NULL AUTO_INCREMENT,
marker varchar(50) NOT NULL,
allele varchar(50),
sequence varchar(500),
avg_no_reads decimal(10,1),
count_seq int,
frequency decimal (5,4),
PRIMARY KEY (id)
);
