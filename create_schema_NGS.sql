
CREATE SCHEMA `NGS_FORENSIC` DEFAULT CHARACTER SET latin1 ;

create table NGS_FORENSIC.Heads(
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



create table NGS_FORENSIC.AutoSTRdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
sequence varchar(500),
CE_validation varchar(50),
head_id int,
PRIMARY KEY (id),
CONSTRAINT head_fk FOREIGN KEY (head_id) REFERENCES Heads(id)
);

CREATE VIEW NGS_FORENSIC.MarkerAutoSTR AS
SELECT m.marker,
       m.allele,
       m.sequence,
       m. avg_no_reads,
       m.count_seq,
       m.count_seq/n.sum_count_seq AS frequency
FROM
  (SELECT DISTINCT (a.sequence), a.allele,
                                 a.marker,
                                 t.count_seq,
                                 t.avg_no_reads
   FROM AutoSTRdata AS a
   JOIN
     (SELECT SEQUENCE,
             marker,
             CE_validation,
             count(SEQUENCE) AS count_seq,
             avg(no_reads) AS avg_no_reads
      FROM AutoSTRdata
      GROUP BY SEQUENCE,
               marker,
               CE_validation) AS t ON a.sequence = t.sequence
   AND a.marker = t.marker
   AND a.CE_validation = t.CE_validation
   WHERE a.CE_validation = 'validated_CE'
     OR a.CE_validation = 'validated_no_reads' ) AS m
JOIN
  (SELECT marker,
          sum(count_seq) AS sum_count_seq
   FROM
     (SELECT DISTINCT (a.sequence), a.allele,
                                    a.marker,
                                    t.count_seq,
                                    t.avg_no_reads
      FROM AutoSTRdata AS a
      JOIN
        (SELECT SEQUENCE,
                marker,
                CE_validation,
                count(SEQUENCE) AS count_seq,
                avg(no_reads) AS avg_no_reads
         FROM AutoSTRdata
         GROUP BY SEQUENCE,
                  marker,
                  CE_validation) AS t ON a.sequence = t.sequence
      AND a.marker = t.marker
      AND a.CE_validation = t.CE_validation
      WHERE a.CE_validation = 'validated_CE'
        OR a.CE_validation = 'validated_no_reads' ) AS x
   GROUP BY marker) AS n ON m.marker = n.marker
ORDER BY marker,
         allele ;
         
         
"or you can create table instead of view and after to do SELECT

create table NGS_FORENSIC.MarkersAutoSTR(
id INT NOT NULL AUTO_INCREMENT,
marker varchar(50) NOT NULL,
allele varchar(50),
sequence varchar(500),
avg_no_reads decimal(10,1),
count_seq int,
frequency decimal (5,4),
PRIMARY KEY (id)
);
