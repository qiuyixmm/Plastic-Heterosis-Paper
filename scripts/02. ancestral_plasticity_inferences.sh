#### Multiple branch estimates
## For ab libitum feeding condition
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/python3.9/envs/CAGEE/bin/cagee --cores 5 \
 --tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/tree.txt \
 --sigma_tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/tree.sigma.alter.txt \
 --infile /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/Pek.Mus.Gos.AL.expression.TPM.txt \
 --output_prefix CAGEE.res.alter

## For overfeeding condition
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/python3.9/envs/CAGEE/bin/cagee --cores 5 \
 --tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/02.OV/tree.txt \
 --sigma_tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/02.OV/tree.sigma.alter.txt \
 --infile /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/02.OV/Pek.Mus.Gos.OV.expression.TPM.txt \
 --output_prefix CAGEE.res.alter
