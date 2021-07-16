package WheatVMap2S1000;

import AoUtils.AoFile;
import AoUtils.CalVCF;
import AoUtils.CountSites;
import WheatGeneticLoad.FilterVCF;
import WheatGeneticLoad.FilterVCF2;
import WheatGeneticLoad.VariantsSum;
import com.google.common.collect.Table;
import daxing.common.IOTool;
import daxing.common.PGF;
import daxing.common.RowTableTool;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.graph.tSaw.TablesawUtils;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;
import tech.tablesaw.api.IntColumn;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.logging.Filter;
import java.util.stream.IntStream;

public class VMap2S1000 {
    public VMap2S1000(){

//        this.bgzip();
//        this.rename();
//        this.vcfQC(); //include many methods XXXXXXX

        /**
         * QC
         */
//        new FilterVCF2().QC();
//        new FilterVCF2().statVcfDepth_SD();
//        new FilterVCF2().mkDepthOfVMapII();
//        new FilterVCF2().mkDepthSummary();
//        this.extractVCF();
        new FilterVCF2().QC();






//        this.geneInfo(); //列出所要建立数据库的基因
//        this.snpAnnotationBuild(); //include many methods


    }

    public void extractVCF(){

        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/008_WheatVMap2_GermplasmInfo_20210708.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subset/chrDsubgenome.100k.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF/diploid.vcf.gz";
        List<String> taxaList = new ArrayList<>();
        RowTable t = new RowTable(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String vcfID = t.getCellAsString(i,0);
            String genomeType = t.getCellAsString(i,4);
            if (genomeType.equals("DD")){
                taxaList.add(vcfID);
            }
        }
        CalVCF.extractVCF(infileS2,outfileS,taxaList);


//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/009_WheatVMap2_GermplasmInfo_20210708.txt";
//        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subset/chrABsubgenome.100k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF/hexaploid.vcf.gz";
//        List<String> taxaList = new ArrayList<>();
//        RowTable t = new RowTable(infileS);
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            String vcfID = t.getCellAsString(i,0);
//            String genomeType = t.getCellAsString(i,4);
//            if (genomeType.equals("AABBDD")){
//                taxaList.add(vcfID);
//            }
//        }
//
//
//        CalVCF.extractVCF(infileS2,outfileS,taxaList);

    }

    public void snpAnnotationBuild(){
//        this.mkGeneVCF();
//        this.mkGeneVCF2();
        this.extractInfoFromGeneVCF();
//        new VariantsSum().extractInfoFromVMap2();

//        new VariantsSum().mkExonAnnotation2();
//        new VariantsSum().addSift();
//        new VariantsSum().addAncestral();
//        new VariantsSum().addDerived_SIFT();
//        new VariantsSum().addDAF();
//        new VariantsSum().addGerp();
//        new VariantsSum().mergeExonSNPAnnotation();
    }


    public void extractInfoFromGeneVCF () {
        /**
         * inputFile
         */
        String outDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/001_genicSNPByChr";
        String vmapDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/103_VMap2.1";
        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/001_geneHC/geneHC.txt";
        /**
         * parameters need to modify
         */


        int subLength = 200;
        File[] fs  = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        AoFile.readheader(geneHCFileS);
        tech.tablesaw.api.Table t = TablesawUtils.readTsv(geneHCFileS);
        System.out.println(t.structure());
        t.sortAscendingOn("Chr", "TranStart");
        IntColumn chrColumn = t.intColumn("chr");
        int chrNum = chrColumn.countUnique(); //chr的个数
        TIntList[] startLists = new TIntList[chrNum]; //list类型的数组，每个数组存放一堆list值
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }

        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3))); //获取开始位点
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4))); //获取终止位点
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1)); //获取基因名字
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_genicSNP.txt")).getAbsolutePath();
            int[] dc = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(dc);
            StringBuilder sb = new StringBuilder();
            if (Arrays.binarySearch(dc, chrIndex+1) < 0) {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB\tTranscript");
            }
            else {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D\tTranscript");
            }
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) { //如果该行的字符长度小于 200，那么最长设置为该行的实际长度
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在开始位点搜索是否在
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    sb.append("\t").append(l.get(4)).append("\t");
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]).append("\t").append(ll.get(7).split("=")[1]).append("\t").append(ll.get(8).split("=")[1]);
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        //在HPC上运行： java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_extractInfoFromVMap2_20200606.txt 2>&1 &
    }


    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkGeneVCF2 () {
        /**
         * inputFile
         */
        String vmapDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF"; //modify
        RowTableTool<String> rowTable = new RowTableTool<>(hcGeneFileS);
        Predicate<List<String>> p= l->Integer.parseInt(l.get(3))==0;
        rowTable.removeIf(p);
        List<String> geneList = rowTable.getColumn(0);
        PGF pgf = new PGF(geneFeatureFileS);
        Predicate<PGF.Gene> predicate = gene -> !geneList.contains(gene);
        pgf.removeIf(predicate);
        pgf.sortGeneByGeneRange();
        List<File> files = IOTool.getFileListInDirEndsWith(vmapDirS, ".vcf.gz");
        String[] outNames=  files.stream().map(File::getName).map(s -> s.replaceAll("_vmap2.0.vcf.gz","_gene_vmap2.0.vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outputDirS, outNames[e]))) {
                String line, subLine;
                int chrID; int pos;
                List<String> temp;
                while ((line= br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    subLine=line.substring(0, 100);
                    temp=PStringUtils.fastSplit(subLine);
                    chrID =Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    int geneIndex=pgf.getGeneIndex(chrID, pos);
                    if (geneIndex < 0) continue;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_mkGeneVCF_20210715.txt 2>&1 &
    }

    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkGeneVCF () {
        /**
         * inputFile
         */
        String vmapDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF"; //modify
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();

        RowTable<String> t = new RowTable<>(hcGeneFileS);
        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        List<Integer> chrList = new ArrayList<>();
        for (int i = 0; i < chrSet.size(); i++) {
            chrList.add(i+1);
        }
        chrList.parallelStream().forEach(chrID -> {
            String inputVCF = new File(vmapDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_vmap2.0.vcf.gz").getAbsolutePath();
            String outputVCF = new File (outputDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_gene_vmap2.0.vcf").getAbsolutePath();
            List<String> geneList = new ArrayList<>();
            List<String> tranList = new ArrayList<>();
            TIntArrayList tranStartList = new TIntArrayList();
            TIntArrayList tranEndList = new TIntArrayList();
            //获取当前 chr 包含的所有 gene 和 trans
            for (int i = 0; i < t.getRowNumber(); i++) {
                int currentChr = Integer.parseInt(t.getCell(i, 5));
                int ifunique = Integer.parseInt(t.getCell(i,3));
                if (ifunique == 0) continue; //过滤不是unique的基因
                if (currentChr < chrID) continue;
                else if (currentChr > chrID) break;
                geneList.add(t.getCell(i, 0));
                tranList.add(t.getCell(i, 4));
                tranStartList.add(t.getCellAsInteger(i,6));
                tranEndList.add(t.getCellAsInteger(i,7));
            }

            //获取每个基因对应的最长转录本的起始终止位点，然后加入大库的 all gene list
            int[] starts = tranStartList.toArray();
            int[] ends = tranEndList.toArray();

            try {
                BufferedReader br = AoFile.readFile(inputVCF);
                BufferedWriter bw = AoFile.writeFile(outputVCF);
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")) {
                    bw.write(temp); bw.newLine();
                }
                bw.write(temp); bw.newLine(); //#CHROM 这一行
                List<String> l = new ArrayList<>();
                int index = -1;
                int pos = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 100));
                    pos = Integer.parseInt(l.get(1));
                    index = Arrays.binarySearch(starts, pos);
                    if (index < 0) index = -index - 2;
                    if (index < 0) continue;
                    if (pos < ends[index]) {
                        cnt++;
                        bw.write(temp);bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(chrID + "\t" + cnt  + "\t" + "mkGeneVCF");

            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_mkGeneVCF_20210715.txt 2>&1 &
    }

    /**
     * 在wheat_v1.1_nonoverlap.txt文件中添加列信息：Chr,TransStart,TransEnd,TranStrand,CDSExonNumber,CDSLength六列信息
     */
    public void geneInfo(){
        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        PGF pgf = new PGF(geneFeatureFileS);
        pgf.sortGeneByName();
        gf.sortGeneByName();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tChr\tTransStart\tTransEnd\tTranStrand\tExonNumber\tCDSLength");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(0);
                String trans = l.get(4);
                int geneIndex = gf.getGeneIndex(gene);
                int longIndex = gf.getLongestTranscriptIndex(geneIndex);
                String longTrans = gf.getTranscriptName(geneIndex,longIndex);
                if (!trans.equals(longTrans)) System.out.println(temp);
                int chr = gf.getChromosomeOfGene(geneIndex);
                if (chr==0)continue;
                int start = gf.getTranscriptStart(geneIndex,longIndex);
                int end = gf.getTranscriptEnd(geneIndex,longIndex);
                int strand = gf.getTranscriptStrand(geneIndex,longIndex);
                int exonNum = gf.getExonList(geneIndex,longIndex).size();
                int CDSlength = pgf.getCDSLen(geneIndex, longIndex);
                sb.append("\t").append(chr).append("\t").append(start).append("\t").append(end).append("\t").append(strand).append("\t").append(exonNum).append("\t").append(CDSlength);
                bw.write(temp + sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }



    public void rename(){
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("mv chr" + chr + ".vmap2.vcf.gz chr" + chr + "_vmap2.0.vcf.gz");
            System.out.println("mv chr" + chr + ".vmap2.vcf.gz.tbi chr" + chr + "_vmap2.0.vcf.gz.tbi");
        }
    }

    /**
     * 对Fei过滤生成的filter2 进行解压缩，用bgzip压缩，再建立索引
     */
    public void bgzip(){
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("nohup gunzip chr" + chr + ".vmap2.vcf.gz 2>&1 &");
            System.out.println("nohup bgzip -@ 4 chr" + chr + ".vmap2.vcf && tabix -p vcf chr" + chr + ".vmap2.vcf.gz &");
        }
    }


}
