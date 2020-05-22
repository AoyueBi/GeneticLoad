package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.SplitScript;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genotype.GenoIOFormat;
import pgl.infra.dna.genotype.GenotypeGrid;
import pgl.infra.dna.genotype.GenotypeOperation;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;

public class FilterVCF2 {

    public FilterVCF2(){

//        this.filter2();
//        this.script();
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_script/sh_filterMAFmissOccurrence20200522.sh",3,11);
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_script/sh_filterMAFmissOccurrence_2_20200522.sh",4,2);

//        this.filter_parallel();

//        String a = "";
//        String b = "";
//        this.filter_singleThread(a,b);

//        this.filterMafbyPopHTD();
//        this.extractField();
        this.getOccurrenceSite();

    }


    /**
     * 检查哪些 MAF 小于 0.01 但是又属于 occurrence 大于 2 （即在2个个体中出现 minor allele 的位点）
      */
    public void getOccurrenceSite(){

        String posMAFfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/chr001_occu2_maf0.01_miss0.2.txt";
        String posOccurrenceFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/occurrence/chr001_occu2_maf0.01_miss0.2.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/001_data/chr001.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/002_occurrenceFile/chr001.vcf";

        TIntArrayList posmaf = AoFile.getTIntList_withoutHeader(posMAFfileS,0);
        TIntArrayList posoccurr = AoFile.getTIntList_withoutHeader(posOccurrenceFileS,0);
        posmaf.sort();
        posoccurr.sort();

        try {

            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(1));
                int index = posmaf.binarySearch(pos);
                int index2 = posoccurr.binarySearch(pos);
                if (index > -1) continue;
                if (index2 < 0) continue;
                bw.write(temp);
                bw.newLine();
                cnt++;
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }


    /**
     * 提取VCF文件的某几列信息
     */
    public void extractField(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_out2_byAomethod";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/004_extractColumn";
        String[] chrArr ={"001","005"};

        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_occu2_maf0.01_miss0.2.vcf").getAbsolutePath();
            String outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();
            int[] column = {1};
            AoFile.extractFileColumn(infileS,"#",column,outfileS);
        }
    }

    public void script(){
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF";
//        String logDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/log";
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};


        String logDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/log/004";

//        String[] chrArr = {"010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        String[] chrArr = {"003","004","005","006","007","008","009"};

        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + ".vcf").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(infileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            System.out.println("java -Xms200g -Xmx200g -jar 044_filterMafOccurrenceMiss.jar " + infileS + " > " + logfileS );
        }
    }

    public void filter_singleThread (String infileS, String outfileDirS) {

//        String outputVCFDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_MAF0.01";
//        String outputVCFDirS = "/data4/home/aoyue/vmap2/analysis/023_hapScanner_basedPopDepth/004_filterMAFmissOccurrence";


        String ABTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/EmmerWheat_S187.txt";
        String ABDTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/BreadWheat_S420.txt";
        String DTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/Ae.tauschii_S36.txt";

        int occu = 2;
        float mafThresh = (float)0.01;
        float missingThresh = (float)0.2;
        File f = new File(infileS);

        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(ABTaxaFileS,0);
        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(ABDTaxaFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(DTaxaFileS,0);

        StringBuilder sb = new StringBuilder();
        sb.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".vcf");
        String outfileS = new File (outfileDirS, sb.toString()).getAbsolutePath();
        int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
        String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);

        //总的VCF文件，获取 gt 和 allTaxa
        GenotypeGrid gt = new GenotypeGrid(f.getAbsolutePath(), GenoIOFormat.VCF);
        gt.sortByTaxa();
        String[] allTaxa = gt.getTaxaNames();
        /**
         * 建立 2 个群体的index
         */
        // 将VCF文件中的两个群体拆开
        int[][] taxaIndices = new int[2][];
        String[][] subTaxa = new String[2][];
        GenotypeGrid[] gts = new GenotypeGrid[2]; //群体的 genotypeTable
        if (subgenome.equals("D")) {
            taxaIndices[0] = new int[abdTaxa.length];
            taxaIndices[1] = new int[dTaxa.length];
            subTaxa[0] = abdTaxa;
            subTaxa[1] = dTaxa;
        }
        else {
            taxaIndices[0] = new int[abTaxa.length];
            taxaIndices[1] = new int[abdTaxa.length];
            subTaxa[0] = abTaxa;
            subTaxa[1] = abdTaxa;
        }

        for (int i = 0; i < taxaIndices.length; i++) { // 第一维，确定 TaxaIndices是哪个群体的
            for (int j = 0; j < taxaIndices[i].length; j++) { //第二维，确定群体内每个 taxa 在VCF文件中的 allTaxa中的索引
                taxaIndices[i][j] = Arrays.binarySearch(allTaxa, subTaxa[i][j]);
                if (taxaIndices[i][j] < 0) System.out.println(subTaxa[i][j]);
            }
        }

        // 从上文建立的 taxaIndex里，选出 单个群体，new gts[i]
        for (int i = 0; i < gts.length; i++) {
            gts[i] = GenotypeOperation.getSubsetGenotypeByTaxon(gt, taxaIndices[i]);
        }
        TIntArrayList posList = new TIntArrayList(); //这是我们将要保留的sites
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            double[] missing = new double[2];
            for (int j = 0; j < 2; j++) { //第一步，获取2个群体的 missing rate
                missing[j] = (double)gts[j].getMissingNumberBySite(i)/gts[j].getTaxaNumber();
            }
            if (missing[0] > missingThresh && missing[1] > missingThresh) continue; //如果缺失率都大于0.2的话，该位点不保留。

            for (int j = 0; j < gts.length; j++) { //判断2个群体的maf值和occurrence 值

                double maf = gts[j].getMinorAlleleFrequency(i);
                if (!(maf < mafThresh)) { //
                    posList.add(gt.getPosition(i));
                    break;
                }
                else if (!(gts[j].getAlternativeAlleleOccurrenceBySite(i) < occu)) {
                    posList.add(gt.getPosition(i));
                    break;
                }
            }
        }
        int[] positions = posList.toArray(); //*********** 该pos库用作下文搜索用   *****************//
        Arrays.sort(positions);
        System.out.println(positions.length + " sites were kept in " + new File(outfileS).getName());
        System.out.println("***********************************************************************");
        System.out.println("***********************************************************************");


        /**
         * begin to write the kept site in vcf file
         */
        try{
            List<Integer> indexABD = new ArrayList<>();
            List<Integer> indexABorD = new ArrayList<>();
            String[] taxaABorDArray = null;
            String aaf = null;  //在注释文件中是写AAF_D 还是AAF_AB
            String annoHeader = null; //Header 中是 四倍体还是六倍体
            if (subgenome.equals("D")){
                taxaABorDArray = dTaxa;
                aaf = "AAF_D";
                annoHeader = this.annotationHeader_Dsub();
            }else{
                taxaABorDArray = abTaxa;
                aaf = "AAF_AB";
                annoHeader = this.annotationHeader_ABsub();
            }

            BufferedReader br = AoFile.readFile(f.getAbsolutePath());
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(annoHeader); bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            while((temp=br.readLine()) != null){
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#CHROM")){
                    l = PStringUtils.fastSplit(temp);
                    bw.write(temp);
                    bw.newLine();
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index1 = Arrays.binarySearch(abdTaxa, taxon);
                        int index2 = Arrays.binarySearch(taxaABorDArray, taxon);
                        if (index1 > -1) {
                            indexABD.add(i);
                        }
                        if (index2 > -1) {
                            indexABorD.add(i);
                        }
                    }
                    Collections.sort(indexABD);
                    Collections.sort(indexABorD);
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    int kk = Arrays.binarySearch(positions,pos);
                    if (kk < 0) continue;
                    String altList = l.get(4);
                    List<String> lgeno = new ArrayList<>();
                    List<String> lABDGeno = new ArrayList<>();
                    List<String> lABorDGeno = new ArrayList<>();

                    for (int i = 9; i < l.size(); i++) {
                        lgeno.add(l.get(i));
                    }
                    for (int i = 0; i < indexABD.size(); i++) {
                        lABDGeno.add(l.get(indexABD.get(i)));
                    }
                    for (int i = 0; i < indexABorD.size(); i++) {
                        lABorDGeno.add(l.get(indexABorD.get(i)));
                    }

                    String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                    String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                    String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                    String INFO = this.getInfo(genoArray, altList);
                    String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
                    String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];


                    StringBuilder sbb = new StringBuilder();
                    for (int i = 0; i < 7; i++) {
                        sbb.append(l.get(i)).append("\t");
                    }

                    sbb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";").append(aaf).append("=").append(ABorDAAF).append("\tGT:AD:GL");
                    for (int i = 9; i < l.size(); i++) {
                        sbb.append("\t").append(l.get(i));
                    }
                    bw.write(sbb.toString());
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(outfileS + " is completed.");

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public String annotationHeader_Dsub(){
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);

        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n" +
                "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
                "##fileDate=" + S.split(" ")[0] + "\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n" +
                "##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n" +
                "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n" +
                "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n" +
                "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n" +
                "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n" +
                "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n" +
                "##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">\n" +
                "##INFO=<ID=AAF_D,Number=1,Type=Float,Description=\"Alternative allele frequency on tetraploid emmer wheat\">\n" +  //这一行根据 D 还是 AB 有变化
                "##ALT=<ID=D,Description=\"Deletion\">\n" +
                "##ALT=<ID=I,Description=\"Insertion\">\n" +
                "##Species=Wheat\n" +
                "##ReferenceGenome=iwgsc_refseqv1.0\n" +
                "##VariantsMapVersion=\"vmap2\"");

        return sb.toString();

    }

    public String annotationHeader_ABsub(){
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);

        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n" +
                "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
                "##fileDate=" + S.split(" ")[0] + "\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n" +
                "##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n" +
                "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n" +
                "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n" +
                "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n" +
                "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n" +
                "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n" +
                "##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">\n" +
                "##INFO=<ID=AAF_AB,Number=1,Type=Float,Description=\"Alternative allele frequency on tetraploid emmer wheat\">\n" +   //这一行根据 D 还是 AB 有变化
                "##ALT=<ID=D,Description=\"Deletion\">\n" +
                "##ALT=<ID=I,Description=\"Insertion\">\n" +
                "##Species=Wheat\n" +
                "##ReferenceGenome=iwgsc_refseqv1.0\n" +
                "##VariantsMapVersion=\"vmap2\"");

        return sb.toString();
    }


    public void filter_parallel () {

//        String inputVCFDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/101_rawMergedVCF";
//        String outputVCFDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/102_MAF0.01";
//        String ABTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/EmmerWheat_S187.txt";
//        String ABDTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/BreadWheat_S420.txt";
//        String DTaxaFileS = "/data4/home/aoyue/vmap2/analysis/000_taxaList/101_taxaListbyPloidy/Ae.tauschii_S36.txt";


        String inputVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/001_data";
//        String outputVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out";
        String outputVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/002_out/occurrence";

        String ABTaxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/EmmerWheat_S187.txt";
        String ABDTaxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/BreadWheat_S420.txt";
        String DTaxaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/010_taxaList/Ae.tauschii_S36.txt";

        int occu = 2;
        float mafThresh = (float)0.01;
        float missingThresh = (float)0.2;
        List<File> fList = AoFile.getFileListInDir(inputVCFDirS);
        Collections.sort(fList);

        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(ABTaxaFileS,0);
        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(ABDTaxaFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(DTaxaFileS,0);

        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".txt");
            String outfileS = new File (outputVCFDirS, sb.toString()).getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);

            //总的VCF文件，获取 gt 和 allTaxa
            GenotypeGrid gt = new GenotypeGrid(f.getAbsolutePath(), GenoIOFormat.VCF_GZ);
            gt.sortByTaxa();
            String[] allTaxa = gt.getTaxaNames();
            /**
             * 建立 2 个群体的index
             */
            // 将VCF文件中的两个群体拆开
            int[][] taxaIndices = new int[2][];
            String[][] subTaxa = new String[2][];
            GenotypeGrid[] gts = new GenotypeGrid[2]; //群体的 genotypeTable
            if (subgenome.equals("D")) {
                taxaIndices[0] = new int[abdTaxa.length];
                taxaIndices[1] = new int[dTaxa.length];
                subTaxa[0] = abdTaxa;
                subTaxa[1] = dTaxa;
            }
            else {
                taxaIndices[0] = new int[abTaxa.length];
                taxaIndices[1] = new int[abdTaxa.length];
                subTaxa[0] = abTaxa;
                subTaxa[1] = abdTaxa;
            }

            for (int i = 0; i < taxaIndices.length; i++) { // 第一维，确定 TaxaIndices是哪个群体的
                for (int j = 0; j < taxaIndices[i].length; j++) { //第二维，确定群体内每个 taxa 在VCF文件中的 allTaxa中的索引
                    taxaIndices[i][j] = Arrays.binarySearch(allTaxa, subTaxa[i][j]);
                    if (taxaIndices[i][j] < 0) System.out.println(subTaxa[i][j]);
                }
            }

            // 从上文建立的 taxaIndex里，选出 单个群体，new gts[i]
            for (int i = 0; i < gts.length; i++) {
                gts[i] = GenotypeOperation.getSubsetGenotypeByTaxon(gt, taxaIndices[i]);
            }
            TIntArrayList posList = new TIntArrayList(); //这是我们将要保留的sites
            for (int i = 0; i < gt.getSiteNumber(); i++) {
                double[] missing = new double[2];
                for (int j = 0; j < 2; j++) { //第一步，获取2个群体的 missing rate
                    missing[j] = (double)gts[j].getMissingNumberBySite(i)/gts[j].getTaxaNumber();
                }
//                if (missing[0] > missingThresh && missing[1] > missingThresh) continue; //如果缺失率都大于0.2的话，该位点不保留。

                for (int j = 0; j < gts.length; j++) { //判断2个群体的maf值和occurrence 值

                    double maf = gts[j].getMinorAlleleFrequency(i);
//                    if (!(maf < mafThresh)) { //
//                        posList.add(gt.getPosition(i));
//                        break;
//                    }
//                    else if (!(gts[j].getAlternativeAlleleOccurrenceBySite(i) < occu)) {
//                        posList.add(gt.getPosition(i));
//                        break;
//                    }

                    if (!(gts[j].getAlternativeAlleleOccurrenceBySite(i) < occu)) {
                        posList.add(gt.getPosition(i));
                        break;
                    }

                }

            }
            int[] positions = posList.toArray();
            Arrays.sort(positions);
            System.out.println(positions.length + " sites were kept in " + new File(outfileS).getName());
            System.out.println("***********************************************************************");
            System.out.println("***********************************************************************");

            /**
             * 测试老师的MAF和我的MAF之间的差别
             *
             */
            try{
                BufferedWriter bw = AoFile.writeFile(outfileS);
                for (int i = 0; i < positions.length; i++) {
                    bw.write(String.valueOf(positions[i]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();

            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }



//            /**
//             * begin to write the kept site in vcf file
//             */
//            try{
//                List<Integer> indexABD = new ArrayList<>();
//                List<Integer> indexABorD = new ArrayList<>();
//                String[] taxaABorDArray = null;
//                String aaf = null;
//                String annoHeader = null;
//                if (subgenome.equals("D")){
//                    taxaABorDArray = dTaxa;
//                    aaf = "AAF_D";
//                    annoHeader = this.annotationHeader_Dsub();
//                }else{
//                    taxaABorDArray = abTaxa;
//                    aaf = "AAF_AB";
//                    annoHeader = this.annotationHeader_ABsub();
//                }
//
//                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
//                BufferedWriter bw = AoFile.writeFile(outfileS);
//                bw.write(annoHeader); bw.newLine();
//                String temp = null;
//                List<String> l = new ArrayList<>();
//                while((temp=br.readLine()) != null){
//                    if (temp.startsWith("##")) continue;
//                    if (temp.startsWith("#CHROM")){
//                        l = PStringUtils.fastSplit(temp);
//                        bw.write(temp);
//                        bw.newLine();
//                        for (int i = 9; i < l.size(); i++) {
//                            String taxon = l.get(i);
//                            int index1 = Arrays.binarySearch(abdTaxa, taxon);
//                            int index2 = Arrays.binarySearch(taxaABorDArray, taxon);
//                            if (index1 > -1) {
//                                indexABD.add(i);
//                            }
//                            if (index2 > -1) {
//                                indexABorD.add(i);
//                            }
//                        }
//                        Collections.sort(indexABD);
//                        Collections.sort(indexABorD);
//                    }
//                    if (!temp.startsWith("#")) {
//                        l = PStringUtils.fastSplit(temp);
//                        int pos = Integer.parseInt(l.get(1));
//                        int kk = Arrays.binarySearch(positions,pos);
//                        if (kk < 0) continue;
//                        String altList = l.get(4);
//                        List<String> lgeno = new ArrayList<>();
//                        List<String> lABDGeno = new ArrayList<>();
//                        List<String> lABorDGeno = new ArrayList<>();
//
//                        for (int i = 9; i < l.size(); i++) {
//                            lgeno.add(l.get(i));
//                        }
//                        for (int i = 0; i < indexABD.size(); i++) {
//                            lABDGeno.add(l.get(indexABD.get(i)));
//                        }
//                        for (int i = 0; i < indexABorD.size(); i++) {
//                            lABorDGeno.add(l.get(indexABorD.get(i)));
//                        }
//
//                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
//                        String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
//                        String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);
//
//                        String INFO = this.getInfo(genoArray, altList);
//                        String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
//                        String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];
//
//
//                        StringBuilder sbb = new StringBuilder();
//                        for (int i = 0; i < 7; i++) {
//                            sbb.append(l.get(i)).append("\t");
//                        }
//
//                        sbb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";").append(aaf).append("=").append(ABorDAAF).append("\tGT:AD:PL");
//                        for (int i = 9; i < l.size(); i++) {
//                            sbb.append("\t").append(l.get(i));
//                        }
//                        bw.write(sbb.toString());
//                        bw.newLine();
//                    }
//                }
//                br.close();
//                bw.flush();
//                bw.close();
//                System.out.println(outfileS + " is completed.");
//
//            }catch (Exception e) {
//                e.printStackTrace();
//                System.exit(1);
//            }
//
        });

    }



//String infileS, String outfileDirS

    public void filterMafbyPopHTD() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/001_data";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/030_FixVMap2/003_out2_byAomethod";

        List<File> fList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fList);

        String hexaFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/BreadWheat_S419.txt";
        String diFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/Ae.tauschii_S36.txt";
        String tetraFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/EmmerWheat_S187.txt";

        int occu = 2;
        float mafThresh = (float) 0.01;
        float missingThresh = (float) 0.2;



        String[] abdTaxa = AoFile.getStringArraybyList_withoutHeader(hexaFileS,0);
        String[] abTaxa = AoFile.getStringArraybyList_withoutHeader(tetraFileS,0);
        String[] dTaxa = AoFile.getStringArraybyList_withoutHeader(diFileS,0);
        Arrays.sort(abdTaxa);
        Arrays.sort(dTaxa);
        Arrays.sort(abTaxa);


//        File f = new File(infileS);

        fList.stream().forEach(f -> {
            StringBuilder s = new StringBuilder();
            s.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".vcf");
            String outfileS = new File(outfileDirS, s.toString()).getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);


            List<Integer> indexABD = new ArrayList<>();
            List<Integer> indexABorD = new ArrayList<>();
            String[] taxaABorDArray = null;
            String aaf = null;  //在注释文件中是写AAF_D 还是AAF_AB
            String annoHeader = null; //Header 中是 四倍体还是六倍体
            if (subgenome.equals("D")) {
                taxaABorDArray = dTaxa;
                aaf = "AAF_D";
                annoHeader = this.annotationHeader_Dsub();
            } else {
                taxaABorDArray = abTaxa;
                aaf = "AAF_AB";
                annoHeader = this.annotationHeader_ABsub();
            }

            System.out.println("Chr\tTotalSNP Num\tBiallelic Num\tTriallelic Num\tDeletion Num\tInsertion Num\tIndel Num");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(annoHeader);
                bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cntSNP = 0; //totalSNP
                int cntkept = 0;

                while ((temp = br.readLine()) != null) {
                    //***********************************************************//
                    if (temp.startsWith("##")) continue;

                    //***********************************************************//

                    if (temp.startsWith("#CHROM")) {
                        l = PStringUtils.fastSplit(temp);
                        bw.write(temp);
                        bw.newLine();

                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index1 = Arrays.binarySearch(abdTaxa, taxon);
                            int index2 = Arrays.binarySearch(taxaABorDArray, taxon);

                            if (index1 > -1) {
                                indexABD.add(i);
                            }
                            if (index2 > -1) {
                                indexABorD.add(i);
                            }
                        }
                        Collections.sort(indexABD);
                        Collections.sort(indexABorD);
                    }
                    if (!temp.startsWith("#")) { //
                        cntSNP++;
                        l = PStringUtils.fastSplit(temp);
                        String altList = l.get(4);
                        int pos  = Integer.parseInt(l.get(1));
                        List<String> lgeno = new ArrayList<>();
                        List<String> lABDGeno = new ArrayList<>();
                        List<String> lABorDGeno = new ArrayList<>();


                        for (int i = 9; i < l.size(); i++) {
                            lgeno.add(l.get(i));
                        }
                        for (int i = 0; i < indexABD.size(); i++) {
                            lABDGeno.add(l.get(indexABD.get(i)));
                        }
                        for (int i = 0; i < indexABorD.size(); i++) {
                            lABorDGeno.add(l.get(indexABorD.get(i)));
                        }

                        String[] genoArray = lgeno.toArray(new String[lgeno.size()]);
                        String[] hexaGenoArray = lABDGeno.toArray(new String[lABDGeno.size()]);
                        String[] ABorDGenoArray = lABorDGeno.toArray(new String[lABorDGeno.size()]);

                        String hexaMAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[1];
                        String diMAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[1];

                        //开始进行,maf判断
                        double hexamaf = Double.parseDouble(hexaMAF);
                        double ABorDmaf = Double.parseDouble(diMAF);

                        if (hexamaf >= 0.01 || (ABorDmaf >= 0.01)) {
                            cntkept++;

                            String INFO = this.getInfo(genoArray, altList);
                            String hexaAAF = this.getSubgenomeInfo(hexaGenoArray, altList).split(",")[0];
                            String ABorDAAF = this.getSubgenomeInfo(ABorDGenoArray, altList).split(",")[0];

                            StringBuilder sb = new StringBuilder();
                            for (int i = 0; i < 7; i++) {
                                sb.append(l.get(i)).append("\t");
                            }
                            sb.append(INFO).append(";AAF_ABD=").append(hexaAAF).append(";").append(aaf).append("=").append(ABorDAAF).append("\tGT:AD:GL");
                            for (int i = 9; i < l.size(); i++) {
                                sb.append("\t").append(l.get(i));
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    } //
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println( cntSNP + "\ttotal\t" + cntkept + " kept is completed at " + outfileS);

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });

    }

    private String getInfo(String[] genoArray, String altList) {
        int dp = 0; //总深度
        int nz = 0; //有基因型的个体数
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的深度统计
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的基因型统计
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型AD还有PL的集合

            //先计算深度
            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
                dp += c; //dp是总深度
                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
            }

            //再计算基因型
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) { //0/0:13,0:0,4,25
                int c = Integer.parseInt(temList.get(j)); // c是基因型0 1 2 其中的一个
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
            int index1 = Integer.parseInt(temList.get(0)); //0/0基因型的
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
            if (index1 != index2) {
                ht++;
            }
        }
        nz = genoArray.length - nz;
        int sum = 0; //所有allele的总数
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }

        StringBuilder sb = new StringBuilder();

        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1); //删除最后一个字符","号
        sb.append(";AC=");
        for (int i = 1; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) { //二维数组的长度是第一维的长度，这里是2（只有1个alt） 或者3 (有2个alt)
            for (int j = i + 1; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";HT=").append(ht).append(";MAF=").append(String.format("%.4f", maf));
        return sb.toString();
    }

    private String getSubgenomeInfo(String[] PopGenoArray, String altList) {
//        int   dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
//        int[] adCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[] acCnt = new int[1 + nAlt]; //所有包括ref和alt的个数
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt]; //GN到底代表什么？
//        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < PopGenoArray.length; i++) {
            if (PopGenoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(PopGenoArray[i], ":"); //tempList是包含基因型AD还有PL的集合
//            temList = PStringUtils.fastSplit(tempList.get(1), ","); //temList是AD所有的深度集合
//            for (int j = 0; j < temList.size(); j++) {
//                int c = Integer.parseInt(temList.get(j)); //c是第j个allele的深度值。注意AD的第一个是ref，第二个是次等位位点的深度，第三个是最小等位位点的深度
//                dp += c; //dp是总深度
//                adCnt[j] += c; //adCnt[j] 是第j个allele的深度值的总和，AD按照 ref alt1 alt2排序
//            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j)); // c是基因型第j个数值
                acCnt[c]++; //acCnt[c] 是所有taxa基因型某一数值如0 1 2的总和
            }
//            int index1 = Integer.parseInt(temList.get(0)); //
//            int index2 = Integer.parseInt(temList.get(1));
//            gnCnt[index1][index2]++; //gnCnt[][]是二维数组，代表alt的个数的矩阵，比如有1个alt，则gnCnt[][]有gnCnt[0][0]  gnCnt[0][1] gnCnt[1][0] gnCnt[1][1]
//            if (index1 != index2) {
//                ht++;
//            }
        }
//        nz = PopGenoArray.length - nz;
        float missRate = (float) ((double) nz/PopGenoArray.length);

        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = (float) ((double) acCnt[1] / sum);
        }
        float aaf = (float) ((double) acCnt[1] / sum);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%.4f", aaf)).append(",").append(String.format("%.4f", maf)).append(",").append(String.format("%.4f",missRate)); //.append(";MAF=")根据实际情况书写MAF_ABD MAF_AB MAF_D
        return sb.toString();
    }
}
