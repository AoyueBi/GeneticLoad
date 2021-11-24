package WheatVMap2S1000;

import AoUtils.AoFile;
import AoUtils.AoMath;
import AoUtils.CalVCF;
import AoUtils.CountSites;
import ExonAnnotation.GeneSNPAnnotation;
import GermplasmInfo.TaxaDB;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class SampleSize2VariantsDiscovery {

    public SampleSize2VariantsDiscovery(){
//        this.mergeExonVCF();
//        this.convert2GenoTable();
//        this.sampleSize2SNPdiscovery();
    }

    /**
     * 估算在 gene 区域，随着样本量的增大，各种类型的变异数目评估，判定是否达到饱和
     */
    public void sampleSize2SNPdiscovery(int loop,int interval, String taxalist,String infileS, String outfileDirS, String exonAnnotationFileS){
//        String[] subspeciesArray = {"ABsub","Dsub"};
//        String[] subspeciesArray = {"Dsub"};

        GeneSNPAnnotation geneSNPAnno = new GeneSNPAnnotation(exonAnnotationFileS);
        System.out.println("******* Stage1: new exonanno class has been build");

        /**
         * 外加一层循环，从取样数1到最大值。 bootstrap
         */
//        int loop = 1;
        for (int j = 0; j < loop; j++) {
            this.mainPipe(j+1,geneSNPAnno,interval,taxalist,infileS,outfileDirS);
            System.out.println("#########*****************************************" + taxalist + "_loop" + loop + "   ####### has been completed.");
//            for (int i = 0; i < subspeciesArray.length; i++) {
//                String choice1 = subspeciesArray[i];
//                this.mainPipe(choice1,j+1,geneSNPAnno);
//                System.out.println("#########*****************************************" + choice1 + "_" + loop + "   ####### has been completed.");
//            }
        }
    }

    /**
     * 根据genotypetable，输出每种分类在每个亚基因组中的个数
     */
    public File getCountinEachCategory(List<String> genoTable,int taxanum,String outfileS,GeneSNPAnnotation geneSNPAnno, int loops){
        //test
//        List<String> l = Arrays.asList("2", "9", "0","2","2","9");
//        String test = this.genotypeCount(l);
//        System.out.println(test);

        //根据类获取分类信息，建立数组
        String[] variantsGroupArray = geneSNPAnno.variantsGroupArray();
        Arrays.sort(variantsGroupArray);
        int[] variantsGroupCount = new int[variantsGroupArray.length];
        //根据是否有ancestral allele,获取derive的 分类信息。建立数组
        String[] loadGroupArray = geneSNPAnno.loadGroupArray();
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
                boolean ifinExonAnnotation = geneSNPAnno.ifinExonAnnotation(id);
                if (!ifinExonAnnotation){ //ifinExonAnnotation 为 true, 不执行；
                    System.out.println(id + " was not in the DB");
                    continue;
                }
                String genotypeCount = this.genotypeCount(geno); //输出形如：0=3；1=2；2=4；9=0 统计各种基因型的个数
                String ifRefisAnc = geneSNPAnno.getIfRefisAnc(id);

                /**
                 * 忽略 ancestral 是否存在，判断和参考基因组0/0相比是否有分离；
                 */
                boolean ifSegration = this.ifSegragation_ignoreIfRefisAnc(genotypeCount);
                if (ifSegration){ //代表有分离，说明检测到了变异的存在
                    String variantsGroup = geneSNPAnno.getVariantsGroup(id);
                    if (variantsGroup == null || variantsGroup == "" || variantsGroup.isEmpty())continue;
                    int index = Arrays.binarySearch(variantsGroupArray,variantsGroup);
                    variantsGroupCount[index]++;
                }

                /**
                 * 如果该位点含有祖先等位基因，并且等于ref或者alt，那么就考虑是否有loadGroup的计数
                 */
                if (!ifRefisAnc.equals("NA")){
                    boolean ifSegragationbyIfRefisAnc = this.ifSegragationbyIfRefisAnc(genotypeCount,ifRefisAnc);
                    if (ifSegragationbyIfRefisAnc){ //代表有分离，说明该位点检测到了derived allele的出现
                        String loadGroup = geneSNPAnno.getLoadGroup(id);
                        if (loadGroup == null || loadGroup == "" || loadGroup.isEmpty())continue;
                        int index = Arrays.binarySearch(loadGroupArray,loadGroup);
//                        if (index<-1)continue;
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

    private boolean ifSegragation_ignoreIfRefisAnc(String genoClassCount){
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

    private boolean ifSegragationbyIfRefisAnc(String genoClassCount,String IfrefisAnc){
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
            if(geno1 >0){ // 0/1
                out = true;
            }else if (geno0>0 && geno2>0){ // 0/0  1/1
                out = true;
            }else if (geno0>0 && geno1>0){ // 0/0 0/1
                out =true;
            }else if (geno1>0 && geno2>0){ //  0/1 1/1
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
    private String genotypeCount (List<String> geno){
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

    public void mainPipe(int loop,GeneSNPAnnotation geneSNPAnno,int interval, String taxalistS,String infileS, String outfileDirS){
        // 非传参法
//        public void mainPipe(String choice1, int loop,GeneSNPAnnotation geneSNPAnno){
//        int interval = Integer.MIN_VALUE; //间隔多少个样品抽一次样
////        String choice = "discreteByConsistentIntervals"; //选择抽样的方式
//        String choice = "discreteByManual";
//        List<String> taxalist = new ArrayList<>();
//        String infileS = null; // genoTable 文件
//        String outfileDirS = null; //结果输出目录
//        String finaloutfileS = null;

        //////// 传参法
        String choice = "discreteByConsistentIntervals"; //选择抽样的方式
//        String choice = "discreteByManual";
        List<String> taxalist = new ArrayList<>();
        taxalist = AoFile.getStringListwithoutHeader(taxalistS,0);

//        if (choice1.equals("Dsub")){ //即抽样6倍体和2倍体
//            interval = 400;
//            taxalist = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/010_sampleSize_to_variantsDiscovery/000_group/TaxainDsub.txt",0);
//            infileS = "/Users/Aoyue/Documents/002_gene_Dsubgenome_genoTable.txt.gz";
//            outfileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/010_sampleSize_to_variantsDiscovery/001/Dsub";
//            finaloutfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/010_sampleSize_to_variantsDiscovery/001/Dsub.txt";
//        }

        //***************************************** 第一阶段： 抽样taxa，非连续抽样 ******************************************************************************************************************************************************//
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
//            int[] sampleArray = {1,100,300,456};
            int[] sampleArray = {1,850};

            List<String>[] taxaListArray =  AoMath.noncontinuousRandom_withoutReplacement(taxaArray,sampleArray);
            taxaListList = Arrays.asList(taxaListArray);
        }


        //********************************************************************************************************************************************************************************************************
        //************************************* 第二阶段：根据第一阶段抽样的taxa,进行genotype table 的提取，并返回数组类型的List,第三阶段：根据geontypeList,进行各种类型的计数 **********************************************************************************************************************************************************//
        List<File> fileList = new ArrayList<>();
        //// 获取基因型table,不输出
        for (int i = 0; i < taxaListList.size(); i++) {
            List<String> taxaList = taxaListList.get(i);
            List<String> genoTable = CalVCF.extractGenotable(infileS,taxaList); //每个table含有header
            int taxanum = taxaList.size(); //当下genotype table含有taxanum个taxa,即当下抽样数
            System.out.println("****** the ID " + (i+1) + " with " + "sample size = " + taxanum + " genotype table has been finished");
            //#### 每个基因型table
            String outfileS = new File(outfileDirS,PStringUtils.getNDigitNumber(3,taxanum) + "_sampleSize_loop" + loop + "_variantsDiscovery.txt").getAbsolutePath();
//            this.getCountinEachCategory(genoTable,taxanum,outfileS,exonanno);
//            fileList.add(this.getCountinEachCategory(genoTable,taxanum,outfileS,geneSNPAnno,loop)); // genoTable是获取的基因型， taxanum是为了写入文件， outfileS是为了写出文件，只有9行，exonanno是个类，为了后续调用，loop也是为了写入文件。
            this.getCountinEachCategory(genoTable,taxanum,outfileS,geneSNPAnno,loop); //不进行文件合并，因为循环多，第四阶段只是合并某一loop 的结果。
            System.out.println("****** the ID " + (i+1) + " with " + "sample size = "+ taxanum + " variants count has been finished");
        }

//        //************************************* 第四阶段：将最终结果进行合并，成为一个文件 **********************************************************************************************************************************************************//
//        File[] fileArray = fileList.toArray(new File[fileList.size()]);
//        String taxa = new File(taxalistS).getName().replaceFirst(".txt","");
//        String finaloutfileS = new File(outfileDirS,taxa+".out.txt").getAbsolutePath(); //### 最终结果放在要输出的文件夹中，并且该文件和 taxalist 的文件名字相同
//        AoFile.mergeTxt_byFileArray(fileArray,finaloutfileS);
//        System.out.println("All done");
    }


    public void convert2GenoTable(){
        String infileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_merge_geneSNPVCF/ABsubgenome.vcf.gz";
        String outfileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/004_geneSNP_genoTable/001_gene_ABsubgenome_genoTable.txt.gz";
        String taxaListFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/004_geneSNP_genoTable/TaxainABsub.txt";
        String[] taxaArray = AoFile.getStringArraybyList_withoutHeader(taxaListFileS,0);
        System.out.println(taxaArray.length);
        CalVCF.extractVCFtable(infileS,taxaArray,outfileS);

        String infileS1 = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_merge_geneSNPVCF/Dsubgenome.vcf.gz";
        String outfileS1 = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/004_geneSNP_genoTable/002_gene_Dsubgenome_genoTable.txt.gz";
        String taxaListFileS1 = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/004_geneSNP_genoTable/TaxainDsub.txt";
        String[] taxaArray1 = AoFile.getStringArraybyList_withoutHeader(taxaListFileS1,0);
        System.out.println(taxaArray1.length);
        CalVCF.extractVCFtable(infileS1,taxaArray1,outfileS1);
    }

    public void mergeExonVCF(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_merge_geneSNPVCF";
        CountSites.mergeVCFtoAB_Dsubgenome(infileDirS,outfileDirS);
    }

}
