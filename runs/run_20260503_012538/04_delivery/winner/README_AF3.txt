AlphaFold 3 抗原-抗体结合测试准备指南

本次筛选出的获胜设计：design_61_0_seq14

该设计包含三个 RSV-F 抗原位点，分别对应以下抗体：

1. 位点 II (Residues 254-277): 对应设计中的第一个 Motif 区域。
   - 推荐抗体：Antibody_Site_II_Palivizumab.fasta (帕利珠单抗)
2. 位点 IV (Residues 422-438): 对应设计中的第二个 Motif 区域。
   - 推荐抗体：Antibody_Site_IV_101F.fasta
3. 位点 V (Residues 163-181): 对应设计中的第三个 Motif 区域。
   - 推荐抗体：Antibody_Site_V_hRSV90.fasta

操作说明：
- 在 AlphaFold 3 提交界面，上传 Winner 目录下的抗原序列。
- 针对每一个位点，分别上传对应抗体重链和轻链序列进行结合面模拟。
- 注意：提供的抗体序列为 Fv 区域（重链和轻链的可变区）。
