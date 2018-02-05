
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

create table NGS_FORENSIC.Heads_Y(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL unique,
project varchar(50) NOT NULL,
analysis varchar(50) NOT NULL,
run varchar(50) NOT NULL,
gender varchar(50) NOT NULL,
created varchar(100) NOT NULL,                             
no_mismatches_Y int,
PRIMARY KEY (id)
);

create table NGS_FORENSIC.Heads_X(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL unique,
project varchar(50) NOT NULL,
analysis varchar(50) NOT NULL,
run varchar(50) NOT NULL,
gender varchar(50) NOT NULL,
created varchar(100) NOT NULL,                             
no_mismatches_X int,
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

create table NGS_FORENSIC.Y_STRdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
sequence varchar(500),
CE_validation varchar(50),
head_Y_id int,
PRIMARY KEY (id),
CONSTRAINT head_Y_fk FOREIGN KEY (head_Y_id) REFERENCES Heads_Y(id)
);

create table NGS_FORENSIC.X_STRdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
sequence varchar(500),
CE_validation varchar(50),
head_X_id int,
PRIMARY KEY (id),
CONSTRAINT head_X_fk FOREIGN KEY (head_X_id) REFERENCES Heads_X(id)
);

create table NGS_FORENSIC.SNPdata(
id INT NOT NULL AUTO_INCREMENT,
sample_name varchar(50) NOT NULL,
marker varchar(50) NOT NULL,
allele varchar(50),
no_reads int,
validation_by_no_reads varchar(50),
head_id int,
PRIMARY KEY (id),
CONSTRAINT head_SNP_fk FOREIGN KEY (head_id) REFERENCES Heads(id)
);

create view NGS_FORENSIC.MarkerAutoSTRview as 
select m.marker, m.allele,  m.sequence,  m. avg_no_reads, m.count_seq, m.count_seq/n.sum_count_seq  as frequency 
from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from AutoSTRdata AS a
			JOIN ( select sequence, marker, CE_validation,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from AutoSTRdata
                group by  sequence, marker, CE_validation) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker and a.CE_validation = t.CE_validation
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as m
join  (select marker,  sum(count_seq) as sum_count_seq
			from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from AutoSTRdata AS a
			JOIN ( select sequence, marker, CE_validation,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from AutoSTRdata
                group by  sequence, marker, CE_validation) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker and a.CE_validation = t.CE_validation
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as x
            group by marker) as n
on m.marker = n.marker 
order by marker, allele ;  

create view NGS_FORENSIC.MarkerY_STRview as 
select m.marker, m.allele,  m.sequence,  m. avg_no_reads, m.count_seq, m.count_seq/n.sum_count_seq  as frequency 
from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from Y_STRdata AS a
			JOIN ( select sequence, marker, CE_validation,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from Y_STRdata
                group by  sequence, marker, CE_validation) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker and a.CE_validation = t.CE_validation
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as m
join  (select marker,  sum(count_seq) as sum_count_seq
			from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from Y_STRdata AS a
			JOIN ( select sequence, marker, CE_validation,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from Y_STRdata
                group by  sequence, marker, CE_validation) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker and a.CE_validation = t.CE_validation
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as x
            group by marker) as n
on m.marker = n.marker 
order by marker, allele ;  

create view NGS_FORENSIC.MarkerX_STRview as 
select m.marker, m.allele,  m.sequence,  m. avg_no_reads, m.count_seq, m.count_seq/n.sum_count_seq  as frequency 
from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from X_STRdata AS a
			JOIN ( select sequence, marker, CE_validation,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from X_STRdata
                group by  sequence, marker, CE_validation) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker and a.CE_validation = t.CE_validation
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as m
join  (select marker,  sum(count_seq) as sum_count_seq
			from (select distinct (a.sequence), a.allele, a.marker, t.count_seq, t.avg_no_reads
			from X_STRdata AS a
			JOIN ( select sequence, marker, CE_validation,
                count(sequence) AS count_seq, avg(no_reads) AS avg_no_reads
                from X_STRdata
                group by  sequence, marker, CE_validation) AS t 
				ON a.sequence = t.sequence AND a.marker = t.marker and a.CE_validation = t.CE_validation
				where  a.CE_validation = 'validated_CE' or a.CE_validation = 'validated_no_reads' ) as x
            group by marker) as n
on m.marker = n.marker 
order by marker, allele ;  

create view NGS_FORENSIC.MarkerSNPview as 
select t1.marker, t1.allele,  t1.avg_no_reads, t1.count_allele, t1.count_allele/t2.sum_count_allele as frequency
from (select t.allele, t.marker, t.validation_by_no_reads, count(t.allele) AS count_allele, avg(t.no_reads) AS avg_no_reads
			from SNPdata as t
			where validation_by_no_reads = 'validated_no_reads'
			group by  allele, marker, validation_by_no_reads
			) as t1
join ( select t1.marker, t1.validation_by_no_reads, sum(t1.count_allele) as sum_count_allele
		from (select t.allele, t.marker, t.validation_by_no_reads, count(t.allele) AS count_allele, avg(t.no_reads) AS avg_no_reads
			from SNPdata as t
			where validation_by_no_reads = 'validated_no_reads'
			group by  allele, marker, validation_by_no_reads
			) as t1
			group by  marker, validation_by_no_reads ) as t2
 on  t1.marker = t2.marker
 order by marker, allele; 




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
