package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.CountSites;
import format.table.RowTable;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TIntArrayList;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class DeleteriousBiologyAoyue {

    public DeleteriousBiologyAoyue() {
//        this.mergeDelgenicSNPAnnotation();
//        this.countDeleteriousVMapII();
//        this.mkDepthOfVMapII();
//        this.mkDepthSummary();
        this.mergeDepthSummary();



    }

    public void mergeDepthSummary(){
        String infileS =  "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Asub.summary.txt";
        String infile2S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Dsub.summary.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
        HashMap<String,String> hm1 = new AoFile().getHashMap(infileS,0,2);
        HashMap<String,String> hm2 = new AoFile().getHashMap(infile2S,0,2);
        Set<String> s = new HashSet<>();
        for (int i = 0; i < hm1.size(); i++) {
            s.addAll(hm1.keySet());
        }
        for (int i = 0; i < hm2.size(); i++) {
            s.addAll(hm2.keySet());
        }
        System.out.println(s.size());

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            String value = null;
            for (int i = 0; i < s.size(); i++) {




            }

            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }




    }

    public void mkDepthSummary () {
        //思想：对每一个taxa做统计，每读进一个taxa，就统计所有深度的和，再除以表格的行数，也即是抽查的位点数，最终得出总得深度除以位点数，得出平均每个位点的深度。
//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Asub";
//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Asub.summary.txt";

        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Dsub";
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_Dsub.summary.txt";


        File[] fs = new File (taxaDepthDirS).listFiles();
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(taxaSummaryFileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                RowTable t = new RowTable (fs[i].getAbsolutePath());
                String taxaName = String.valueOf(t.getHeader().get(0)).replaceFirst("_siteDepth", "");
                double value = 0;
                for (int j = 0; j < t.getRowNumber(); j++) {
                    value+=t.getCellAsDouble(j, 0);
                }
                double dd = (double)value/t.getRowNumber();
                bw.write(taxaName+"\t"+String.valueOf(i+1)+"\t"+String.format("%.2f", dd));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkDepthOfVMapII(){
//        String vcfFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/024_subsetVCF_maf0.01byPop/002_merged/chr.Asub.maf0.01byPop.vcf.gz";
//        String hmpInfoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_A.SNP_anno.txt.gz";
//        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Asub";

        String vcfFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/024_subsetVCF_maf0.01byPop/002_merged/chr.Dsub.maf0.01byPop.vcf.gz";
        String hmpInfoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/014_merge/chr_D.SNP_anno.txt.gz";
        String taxaDepthDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxa_Dsub";

        int snpNum = 0;
        int size = 500000;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(hmpInfoFileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
            }
            snpNum = cnt;
            int[] indices = new int[size];
            for (int i = 0; i < size; i++) {
                indices[i] = (int)(Math.random()*snpNum);
            }
            Arrays.sort(indices);
            br = IOUtils.getTextGzipReader(vcfFileS);
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l = PStringUtils.fastSplit(temp, "\t");
            String[] taxa = new String[l.size()-9];
            for (int i = 0; i < taxa.length; i++) {
                taxa[i] = l.get(i+9);
            }
            TIntArrayList[] depthList = new TIntArrayList[taxa.length];
            for (int i = 0; i < taxa.length; i++) depthList[i] = new TIntArrayList();
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++; // 对snp开始计数
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines");
                int idx = Arrays.binarySearch(indices, cnt-1);
                if (idx < 0) continue;
                l = PStringUtils.fastSplit(temp, "\t");
                for (int i = 0; i < taxa.length; i++) {
                    String genoS = l.get(i+9);
                    if (genoS.startsWith(".")) {
                        depthList[i].add(0);
                        continue;
                    }
                    List<String> ll = PStringUtils.fastSplit(genoS, ":");
                    List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                    int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1));
                    depthList[i].add(depth);
                }
            }
            for (int i = 0; i < taxa.length; i++) {
                String outfileS = new File (taxaDepthDirS, PStringUtils.getNDigitNumber(3, i+1)+"depth.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(taxa[i]+"_siteDepth");
                bw.newLine();
                int[] depth = depthList[i].toArray();
                for (int j = 0; j < depth.length; j++) {
                    bw.write(String.valueOf(depth[j]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void countDeleteriousVMapII() {
        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/del";
        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz";
        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveDeleterious_vmap2.txt";
        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/reccesiveDeleterious_vmap2.txt";
        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }

        RowTable<String> t = new RowTable(deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        new AoFile().readheader(deleFileS);
        String derivedAllele = null;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引
            String ancestralAllele = t.getCell(i, 14);
            String majorAllele = t.getCell(i, 5);
            String minorAllele = t.getCell(i, 6);
            if (ancestralAllele.equals(majorAllele)) {
                derivedAllele = minorAllele;
                posList[index].add(t.getCellAsInteger(i, 2));
                charList[index].add(derivedAllele.charAt(0));
            }
            if (ancestralAllele.equals(minorAllele)) {
                derivedAllele = majorAllele;
                posList[index].add(t.getCellAsInteger(i, 2));
                charList[index].add(derivedAllele.charAt(0));
            }
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
            Arrays.sort(delePos[i]);
        }

        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/004_vmap2.1_TaxaList/TaxaList.txt";
        List<String> taxaList = new AoFile().getStringListwithoutHeader(vmap2TaxaList, 0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);

        int taxaNum = taxa.length;
        double[] addCount = new double[taxa.length];
        int[] recCount = new int[taxa.length];
        int[] siteWithMinDepthCount = new int[taxa.length];

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_del_vmap2.1.vcf.gz";
            delVmapFileS = new File(delVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {

                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
                        }
                        for (int i = 0; i < taxainVCFfile.length; i++) {
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }
                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 1000 == 0) {
                            System.out.println(String.valueOf(cnt) + " lines on chr " + String.valueOf(chr));
                        }

                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1));
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        }
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0;
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                        }
                        idx[1] = 0;

                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            if (v1 == idx[0]) {
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                addCount[taxaIndex] += 0.5;
                            } else {
                                addCount[taxaIndex] += 1;
                                recCount[taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[taxaIndex]++;
                        }

                    }

                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            try {
                BufferedWriter bw = IOUtils.getTextWriter(addCountFileS);
                bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
                bw.newLine();
                for (int i = 0; i < addCount.length; i++) {
                    bw.write(taxa[i] + "\t" + String.valueOf(addCount[i]) + "\t" + String.valueOf(siteWithMinDepthCount[i]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                bw = IOUtils.getTextWriter(recCountFileS);
                bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth");
                bw.newLine();
                for (int i = 0; i < recCount.length; i++) {
                    bw.write(taxa[i] + "\t" + String.valueOf(recCount[i]) + "\t" + String.valueOf(siteWithMinDepthCount[i]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
    }

    public void mergeDelgenicSNPAnnotation() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz";
        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");
        //Total lines without header count is 91306 at merged file /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/del_merged/chr_SNP_anno_del.txt.gz
        //12/26/19 8:37 PM
    }

}
