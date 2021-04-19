package WheatGeneticLoad;

import AoUtils.*;
import ExonAnnotation.ExonAnnotation;
import GermplasmInfo.TaxaDB;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;

public class VMap2Cal {

    public VMap2Cal(){
//        this.extractAAFfromVMap2(); //计算在一个vcf中，不同给定群体的AAF值是多少， linux 完成
//        this.buildJointSFS();
//        this.sampleSize2variantsDiscovery();
//        File out = WheatUtils.getCentromereFile("/Users/Aoyue/Documents/test.txt");

    }


    /**
     * 估算在exon区域，随着样本量的增大，各种类型的变异数目评估，判定是否达到饱和
     */
    public void sampleSize2variantsDiscovery(){
//        this.mergeExonVCF();
//        this.convert2GenoTable();
//        this.mainPipe();

//        String[] subspeciesArray = {"Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid","Landrace-AB","Cultivar-AB","Landrace-D","Cultivar-D","OtherTetraploid","OtherHexaploid-AB","OtherHexaploid-D"};
//        String[] subspeciesArray = {"Landrace","Cultivar","OtherHexaploid"};

        String[] subspeciesArray = {"ABsub","Dsub"};

        /**
         * 外加一层循环，从取样数1到最大值。 每个取样数重复多次
         */

        int loop = 10;
//        for (int j = 0; j < loop; j++) {
//            for (int i = 0; i < subspeciesArray.length; i++) {
//                String choice1 = subspeciesArray[i];
//                this.mainPipe(choice1,j+1);
//                System.out.println("#########*****************************************" + choice1 + "_" + loop + "   ####### has been completed.");
//            }
//        }

//        List<File> fileList = AoFile.getFileListInDir("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/006_sampleSize2VariantsDiscovery_ABsub");
//        File[] fileArray = fileList.toArray(new File[fileList.size()]);
//        AoFile.mergeTxt_byFileArray(fileArray,"/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/007_merge006/ABsub.txt");
//        System.out.println("All done");

        List<File> fileList = AoFile.getFileListInDir("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/006_sampleSize2VariantsDiscovery_Dsub");
        File[] fileArray = fileList.toArray(new File[fileList.size()]);
        AoFile.mergeTxt_byFileArray(fileArray,"/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/007_merge006/Dsub.txt");
        System.out.println("All done");


        //Ae.tauschii, Cultivar, Domesticated_emmer, Free_threshing_tetraploid, Landrace, OtherHexaploid, OtherTetraploid, Wild_emmer
    }

    /**
     * 根据genotypetable，输出每种分类在每个亚基因组中的个数
     */
    public File getCountinEachCategory(List<String> genoTable,int taxanum,String outfileS,ExonAnnotation exonanno, int loops){
        //test
//        List<String> l = Arrays.asList("2", "9", "0","2","2","9");
//        String test = this.genotypeCount(l);
//        System.out.println(test);

        //根据类获取分类信息，建立数组
        String[] variantsGroupArray = exonanno.variantsGroupArray();
        Arrays.sort(variantsGroupArray);
        int[] variantsGroupCount = new int[variantsGroupArray.length];
        //根据是否有ancestral allele,获取derive的 分类信息。建立数组
        String[] loadGroupArray = exonanno.loadGroupArray();
        Arrays.sort(loadGroupArray);
        int[] loadGroupCount = new int[loadGroupArray.length];

        try {
            //根据genoTable 获取 geno 和 chr pos
            String header = genoTable.get(0);
            List<String> l = new ArrayList<>();
            //对每一行进行一个处理
            for (int i = 1; i < genoTable.size(); i++) { //对每行位点进行处理
                String temp = genoTable.get(i); //第一行的内容
                l = PStringUtils.fastSplit(temp);
                String chr = l.get(0);
                String pos = l.get(1);
                List<String> geno = new ArrayList<>();
                for (int j = 2; j < l.size(); j++) { //从第3个元素开始有基因型输出 index2
                    geno.add(l.get(j));
                }
                String id = chr + "-" + pos;
                boolean ifinExonAnnotation = exonanno.ifinExonAnnotation(id);
                if (!ifinExonAnnotation){ //ifinExonAnnotation 为 true, 不执行；
                    continue;
                }
                String genotypeCount = this.genotypeCount(geno); //输出形如：0=3；1=2；2=4；9=0 统计各种基因型的个数
                String ifRefisAnc = exonanno.getIfRefisAnc(id);

                /**
                 * 忽略 ancestral 是否存在，判断和参考基因组0/0相比是否有分离；
                 */
                boolean ifSegration = this.ifSegragation_ignoreIfRefisAnc(genotypeCount);
                if (ifSegration){ //代表有分离，说明检测到了变异的存在
                    String variantsGroup = exonanno.getVariantsGroup(id);
                    int index = Arrays.binarySearch(variantsGroupArray,variantsGroup);
                    variantsGroupCount[index]++;
                }

                /**
                 * 如果该位点含有祖先等位基因，并且等于ref或者alt，那么就考虑是否有loadGroup的计数
                 */
                if (!ifRefisAnc.equals("NA")){
                    boolean ifSegragationbyIfRefisAnc = this.ifSegragationbyIfRefisAnc(genotypeCount,ifRefisAnc);
                    if (ifSegragationbyIfRefisAnc){ //代表有分离，说明该位点检测到了derived allele的出现
                        String loadGroup = exonanno.getLoadGroup(id);
                        int index = Arrays.binarySearch(loadGroupArray,loadGroup);
                        loadGroupCount[index]++;
                    }
                }
            } //所有位点处理完毕

            //开始写文件
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write( "TaxaNum\tGroup\tObservedCount\tLoop");bw.newLine();

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < variantsGroupArray.length; i++) {
                sb.setLength(0);
                sb.append(String.valueOf(taxanum)).append("\t").append(variantsGroupArray[i]).append("\t").append(variantsGroupCount[i]).append("\t").append(loops);
                bw.write(sb.toString());
                bw.newLine();
            }

            for (int i = 0; i < loadGroupArray.length; i++) {
                sb.setLength(0);
                sb.append(String.valueOf(taxanum)).append("\t").append(loadGroupArray[i]).append("\t").append(loadGroupCount[i]).append("\t").append(loops);
                bw.write(sb.toString());
                bw.newLine();
            }

            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        File out = new File(outfileS);
        return out;
    }



    public boolean ifSegragation_ignoreIfRefisAnc(String genoClassCount){
        boolean out = false;
        List<String> l = PStringUtils.fastSplit(genoClassCount,";");
        int geno0 = Integer.parseInt(l.get(0).split("=")[1]);
        int geno1 = Integer.parseInt(l.get(1).split("=")[1]);
        int geno2 = Integer.parseInt(l.get(2).split("=")[1]);
        int geno9 = Integer.parseInt(l.get(3).split("=")[1]);
        if (geno1 >0 || geno2 >0){ //都是1/1，或者有 0/1 代表有分离
            out = true;
        }else out = false;

        return out;
    }

    public boolean ifSegragationbyIfRefisAnc(String genoClassCount,String IfrefisAnc){
        boolean out = false;
        List<String> l = PStringUtils.fastSplit(genoClassCount,";");
        int geno0 = Integer.parseInt(l.get(0).split("=")[1]);
        int geno1 = Integer.parseInt(l.get(1).split("=")[1]);
        int geno2 = Integer.parseInt(l.get(2).split("=")[1]);
        int geno9 = Integer.parseInt(l.get(3).split("=")[1]);
        if (IfrefisAnc.equals("Anc")){ //1/1 或者 0/1 代表有分离
            if (geno1 >0 || geno2 >0){
                out = true;
            }else out = false;
        }
        if (IfrefisAnc.equals("Der")){ // 0/0 0/1 代表有分离
//            if (geno0 >0 || geno1 >0){
//                out = true;
//            }else out = false;

            /**
             * 程序修改，满足该位点必须有分离，并且含有derived 才算是可以发现变异;分为4种情况
             */
            if(geno1 >0){
                out = true;
            }else if (geno0>0 && geno2>0){
                out = true;
            }else if (geno0>0 && geno1>0){
                out =true;
            }else if (geno1>0 && geno2>0){
                out=true;
            }
        }

        return out;
    }

    /**
     * 判断某个位点中，多个个体的基因型中的 0/0 1/1 0/1 ./.个数
     * 对应码：
     * 0/0 -> 0
     * 0/1 -> 1
     * 1/1 -> 2
     * ./. -> 9
     * @return
     */
    public String genotypeCount (List<String> geno){
        String out = null;
        String[] genoClass = {"0","1","2","9"};
        int[] genoClassCount = new int[genoClass.length];
        Arrays.sort(genoClass);
        for (int i = 0; i < geno.size(); i++) {
            String query = geno.get(i);
            int index = Arrays.binarySearch(genoClass,query);
            if (index <0) System.out.println("There is some error about the genotype, please check your genotype.");
            genoClassCount[index]++;
        }

        StringBuilder sb = new StringBuilder();
        sb.append("0=").append(genoClassCount[0]).append(";").append("1=").append(genoClassCount[1]).append(";").
                append("2=").append(genoClassCount[2]).append(";").append("9=").append(genoClassCount[3]);
        out = sb.toString();
        return out;
    }

    public void mainPipe(String choice1, int loop){

        String parentDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/004_sampleSize_to_variantsDiscovery";
        String taxaListFileS = null;
        String infileS = null;
        String outfileDirS = null;
        String outfileDirS2 = null;
        String finaloutfileS = null;
        int interval = Integer.MIN_VALUE;
        List<String> taxalist = new ArrayList<>();
        TaxaDB taxadb = new TaxaDB();



        String[] goal = {"Dsub","ABsub"}; //因为不同亚基因组的样本数不同，所以要分开计算
//        String choice1 = "Ae.tauschii";

        if (choice1.equals("Dsub")){ //即抽样6倍体和2倍体
            interval = 50;
            taxalist = taxadb.getTaxaListbyDsub();
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
//            outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/002_sampleSize2VariantsDiscovery";
            outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/006_sampleSize2VariantsDiscovery_Dsub";
            outfileDirS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/006_sampleSize2VariantsDiscovery_Dsub";
//            finaloutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/003_merge/test.txt";
            finaloutfileS = new File(outfileDirS2,choice1 + "_loop" + loop + ".txt").getAbsolutePath();
        }

        if (choice1.equals("ABsub")){ //即抽样6倍体和4倍体
            interval = 50;
            taxalist = taxadb.getTaxaListbyABsub();
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
//            outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/002_sampleSize2VariantsDiscovery_AB";
            outfileDirS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/006_sampleSize2VariantsDiscovery_ABsub";
//            finaloutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/003_merge/test.txt";
            finaloutfileS = new File(outfileDirS2,choice1 + "_loop" + loop + ".txt").getAbsolutePath();
        }


        if (choice1.equals("Ae.tauschii")){ //即抽样二倍体
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

        if (choice1.equals("Wild_emmer")){ //即抽样四倍体
            interval =5;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

        if (choice1.equals("Domesticated_emmer")){ //即抽样四倍体
            interval =5;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

        if (choice1.equals("Free_threshing_tetraploid")){ //即抽样四倍体
            interval =5;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

        if (choice1.equals("OtherTetraploid")){ //即抽样四倍体
            interval =5;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

//        if (choice1.equals("Landrace")){ //即抽样六倍体
//            interval =5;
//            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
//            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
//            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
//            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
//        }
//
//        if (choice1.equals("Cultivar")){ //即抽样六倍体
//            interval = 2;
//            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
//            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
//            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
//            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
//        }
//
//        if (choice1.equals("OtherHexaploid")){ //即抽样六倍体
//            interval =5;
//            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
//            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
//            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
//            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
//        }

        if (choice1.equals("Landrace")){ //即抽样六倍体
            interval =5;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

        if (choice1.equals("Cultivar")){ //即抽样六倍体
            interval =2;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }

        if (choice1.equals("OtherHexaploid")){ //即抽样六倍体
            interval =5;
            taxalist = taxadb.getTaxaArrayofSubspecies(choice1);
            infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
            outfileDirS = new File(parentDirS,choice1).getAbsolutePath(); new File(outfileDirS).mkdirs();
            finaloutfileS = new File(parentDirS,choice1 + ".txt").getAbsolutePath();
        }


        String[] sampleType = {"continuous","discreteByConsistentIntervals","discreteByManual"};
        String choice = "discreteByConsistentIntervals";

        //第一阶段： 抽样taxa，非连续抽样
//        String[] taxaArray = AoFile.getStringArraybyList(taxaListFileS,0); //本行适用于 AB sub  D sub, 不适用于亚群
        String[] taxaArray = taxalist.toArray(new String[taxalist.size()]);
        List<List<String>> taxaListList = new ArrayList<>();

        //##**************************** 连续抽样
        if (choice.equals("continuous")){
            List<String>[] taxaListArray =  AoMath.continuousRandom_withoutReplacement(taxaArray,10);
            taxaListList = Arrays.asList(taxaListArray);
        }

        //##**************************** 非连续抽样，程序化间隔，即间隔式抽样
        if (choice.equals("discreteByConsistentIntervals")){
            TIntArrayList sampleCountList = new TIntArrayList();
//            int interval = 10;
            sampleCountList.add(1);
            for (int i = 1; i < taxaArray.length; i++) { //每隔100个抽一次样******** 合计抽了6次 1，100，200，300，400，456
                int element = interval*i;
                if (element < taxaArray.length){
                    sampleCountList.add(element);
                }else {
                    sampleCountList.add(taxaArray.length);
                    break;
                }
            }
            List<String>[] taxaListArray =  AoMath.noncontinuousRandom_withoutReplacement(taxaArray,sampleCountList);
            taxaListList = Arrays.asList(taxaListArray);

        }

        //##**************************** 非连续抽样，自定义抽样，即间隔式抽样
        if (choice.equals("discreteByManual")){
            int[] sampleArray = {1,100,300,456};
            List<String>[] taxaListArray =  AoMath.noncontinuousRandom_withoutReplacement(taxaArray,sampleArray);
            taxaListList = Arrays.asList(taxaListArray);
        }


        ////第二阶段：根据第一阶段抽样的taxa,进行genotype table 的提取，并返回数组类型的List,第三阶段：根据geontypeList,进行各种类型的计数
        List<File> fileList = new ArrayList<>();
        //// 获取基因型table,不输出
        ExonAnnotation exonanno = new ExonAnnotation();
        System.out.println("new exonanno class has been build");
        for (int i = 0; i < taxaListList.size(); i++) {
            List<String> taxaList = taxaListList.get(i);
            List<String> genoTable = CalVCF.extractGenotable(infileS,taxaList); //每个table含有header
            System.out.println("the " + String.valueOf(i+1) + " genotype table has been finished");
            int taxanum = taxaList.size(); //当下genotype table含有taxanum个taxa,即当下抽样数
            //#### 每个基因型table
            String outfileS = new File(outfileDirS,PStringUtils.getNDigitNumber(3,taxanum) + "_sampleSize_loop" + loop + "_variantsDiscovery.txt").getAbsolutePath();
//            this.getCountinEachCategory(genoTable,taxanum,outfileS,exonanno);
            fileList.add(this.getCountinEachCategory(genoTable,taxanum,outfileS,exonanno,loop)); // genoTable是获取的基因型， taxanum是为了写入文件， outfileS是为了写出文件，只有9行，exonanno是个类，为了后续调用，loop也是为了写入文件。
            System.out.println("the " + String.valueOf(i+1) + " variants count has been finished");
        }

//        //////第四阶段：将最终结果进行合并，成为一个文件
        File[] fileArray = fileList.toArray(new File[fileList.size()]);
        AoFile.mergeTxt_byFileArray(fileArray,finaloutfileS);
        System.out.println("All done");

        //test
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/test.txt";
//        String[] taxaArray = {"BaiMaZha","BaiQiuMai"};
//        CalVCF.extractGenotable(infileS,taxaArray,outfileS);

        //test
//        List<String> taxaList = new ArrayList<>(); taxaList.add("BaiMaZha");taxaList.add("BaiQiuMai");
//        CalVCF.extractGenotable(infileS,taxaList,outfileS);

    }


    /**
     * 将exonVCF 转化为 derived allele table
     * 假设 ancestral allele 是 0/0 -> 0, derived allele 是1/1 -> 2, heter 是 0/1 ->2
     * coode: 0: ancestral allele 1:derived allele 2:derived
     */
    public void convert2GenoTable(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/ABsubgenome.vcf";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_ABsubgenome_genoTable.txt.gz";
//        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/000_ABsub.txt";
//        String[] taxaArray = AoFile.getStringArraybyList(taxaListFileS,0);
//        System.out.println(taxaArray.length);
//        CalVCF.extractVCFtable(infileS,taxaArray,outfileS);


        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/Dsubgenome.vcf";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exon_genoTable/001_exon_Dsubgenome_genoTable.txt.gz";
        String taxaListFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/000_Dsub.txt";
        String[] taxaArray = AoFile.getStringArraybyList(taxaListFileS,0);
        System.out.println(taxaArray.length);
        CalVCF.extractVCFtable(infileS,taxaArray,outfileS);

    }

    public void mergeExonVCF(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/001_exon.vcf";
//        CalVCF.mergeVCF(infileDirS,outfileS);


        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonVCF/";

        CountSites.mergeVCFtoAB_Dsubgenome(infileDirS,outfileDirS);

    }

    /**
     * 根据R 导入的数据，进行矩阵建立，最终画出2D-SFS图形
     */
    public void buildJointSFS(){

        int[][] matrix = new int[420*2][187*2];
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/008_dafCount.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/rscript/referenceEvaluation/data/009_2dSFS.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            System.out.println(matrix.length);
            System.out.println(matrix[0].length);
            System.out.println();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                int abd = Integer.parseInt(l.get(0));
                int ab = Integer.parseInt(l.get(1));
                matrix[abd-1][ab-1]++;
            }
            br.close();


            bw.write("\t");
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < matrix[0].length; i++) {
                sb.append("AB").append(i+1).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();

            for (int i = 0; i < matrix.length; i++) {
                bw.write("ABD" + (i+1));
                for (int j = 0; j < matrix[0].length; j++) {
                    int count = matrix[i][j];
                    bw.write("\t" + count);
                }
                bw.newLine();
            }

            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * Goal: 从VMap2.1 文件中提取第 8 列 （index 7），并判断maf值的大小
     */
    public void extractAAFfromVMap2(){
        //zcat chr1A_vmap2.1_hexaploid.vcf.gz |head -n 100 |grep -v "#" | awk '{print $8}'|cut -d";" -f7|cut -d"=" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l

//                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};

//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

//        int[] darray = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println( "zcat chr" + chr + "_vmap2.1_hexaploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l > ../test/" + chr + "_lessthan0.05.txt 2>&1 &");
//            System.out.println( "zcat chr" + chr + "_vmap2.1_hexaploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 > 0.05) print $0}'|wc -l > ../test2/" + chr + "_morethan0.05.txt 2>&1 &");

//            System.out.println( "zcat chr" + chr + "_vmap2.1_tetraploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l > ../test3/" + chr + "_lessthan0.05.txt 2>&1 &");
//            System.out.println( "zcat chr" + chr + "_vmap2.1_tetraploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 > 0.05) print $0}'|wc -l > ../test4/" + chr + "_morethan0.05.txt 2>&1 &");


//            System.out.println( "zcat chr" + chr + "_vmap2.1_diploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 <= 0.05) print $0}'|wc -l > ../test5/" + chr + "_lessthan0.05.txt 2>&1 &");
            System.out.println( "zcat chr" + chr + "_vmap2.1_diploid.vcf.gz|grep -v \"#\" | awk '{print $8}'|cut -d\";\" -f7|cut -d\"=\" -f2| awk '{if ($0 > 0.05) print $0}'|wc -l > ../test6/" + chr + "_morethan0.05.txt 2>&1 &");

        }
    }

}
