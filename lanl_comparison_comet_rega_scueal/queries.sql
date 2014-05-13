--
-- total number of pol sequence in the LANL dataset
--
SELECT count(id) FROM sequences WHERE gene = 'POL' AND LENGTH(sequence) > 800;

--
-- referenced subtypes in LANL pol dataset 
-- n = 106634
--
SELECT subtype_lanl,count(id) as anzahl 
FROM sequences 
WHERE gene = 'POL' 
AND LENGTH(sequence) > 800
AND (subtype_LANL LIKE '%\_%' OR subtype_LANL IN ('A','A1','B','C','D','F','F1','F2','G','01_AE','02_AG','H','J','K ','A2'))
GROUP BY subtype_lanl
ORDER BY anzahl DESC;

SELECT count(id)
FROM sequences 
WHERE gene = 'POL' 
AND LENGTH(sequence) > 800
AND (subtype_LANL LIKE '%\_%' OR subtype_LANL IN ('A','A1','B','C','D','F','F1','F2','G','01_AE','02_AG','H','J','K ','A2'));

--
-- all referenced subtype apart 'A','A1','B','C','D','F','F1','F2','G','01_AE','02_AG'
-- n = 3530
--
SELECT count(id)
FROM sequences 
WHERE gene = 'POL' 
AND LENGTH(sequence) > 800
AND (subtype_LANL LIKE '%\_%' OR subtype_LANL IN ('A','A1','B','C','D','F','F1','F2','G','01_AE','02_AG','H','J','K ','A2'))
AND subtype_LANL NOT IN ('A','A1','B','C','D','F','F1','F2','G','01_AE','02_AG');

