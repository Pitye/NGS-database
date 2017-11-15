
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


