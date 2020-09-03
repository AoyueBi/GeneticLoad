package AoUtils;

/**
 *
 */

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * @author AoyueBi
 *
 */
public class AoString {
    public AoString(){

    }

    /**
     * 根据输入文件夹名称，自动构建输出文件夹名称，
     * 即将/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/002_0.0001_100_500
     * 改成/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/102_0.0001_100_500
     * @param infileDirS
     * @return
     */
    public static String autoOutfileDirS (String infileDirS){
        String out = null;
        String inDirSName = new File(infileDirS).getName();
        char[] cc = new char[inDirSName.length()];

        char first = inDirSName.charAt(0);
        int a = Character.getNumericValue(first);
        char first_modify = String.valueOf(a+1).charAt(0);

        cc[0] = first_modify;
        for (int i = 1; i < cc.length; i++) {
            cc[i] = inDirSName.charAt(i);
        }
        out = new String(cc);

        out = new File(new File(infileDirS).getParent(),out).getAbsolutePath();
        return out;
    }

    /**
     * 将 double类型的数据格式化输出，保留4位有效数字
     * @param value
     * @return
     */
    public static String outFormatted(double value){
        String out = null;
        out = String.format("%.4f",value);

        return out;
    }

    /**
     * 根据 001 002 or 005 判断该染色体是属于 AB四倍体还是D二倍体
     * @param chr
     * @return
     */
    public static String ABorD(String chr){
        String out = null;
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
        Arrays.sort(chrArr);
        int index = Arrays.binarySearch(chrArr,chr);
        if (index > -1) out = "d";
        if (index < 0) out = "ab";
        return out;
    }


    /**
     * 根据基因的名字获取染色体信息 //TraesCS1D02G011700
     * @param gene
     * @return
     */
    public static String getChrFromGene(String gene){ return gene.substring(7,9);}

    /**
     * 根据chr名字获取亚基因组信息 //001
     * @param chr
     * @return
     */
    public static String getSubfromChrID(int chr){
        String out = null;
//        String[] chrAArr ={"001","002","007","008","013","014","019","020","025","026","031","032","037","038"};
//        String[] chrBArr ={"003","004","009","010","015","016","021","022","027","028","033","034","039","040"};
//        String[] chrDArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        int[] chrAArr ={1,2,7,8,13,14,19,20,25,26,31,32,37,38};
        int[] chrBArr ={3,4,9,10,15,16,21,22,27,28,33,34,39,40};
        int[] chrDArr ={5,6,11,12,17,18,23,24,29,30,35,36,41,42};

        Arrays.sort(chrAArr); Arrays.sort(chrBArr); Arrays.sort(chrDArr);

        int indexA = Arrays.binarySearch(chrAArr,chr);
        int indexB = Arrays.binarySearch(chrBArr,chr);
        int indexD = Arrays.binarySearch(chrDArr,chr);

        if (indexA > -1) out = "A";
        if (indexB > -1) out = "B";
        if (indexD > -1) out = "D";
        return out;
    }

    /**
     * 根据转录本的名字获取亚基因组信息 //TraesCS1A02G001800.1
     * @param gene
     * @return
     */
    public static String getSubFromTranscript(String gene){ return gene.substring(8,9);}

    /**
     * 主要用于脚本书写，在hapScanner 脚本中用到
     * @param chrArr
     * @return
     */
    public static String getPloidy(String[] chrArr){
        String out = null;
        if (chrArr.length == 42) out = "abd";
        if (chrArr.length == 28) out = "ab";
        if (chrArr.length == 14) out = "d";
        return  out;
    }
//String name = new File(outPath).getName();
    //String chr = name.substring(name.indexOf("chr")+3,name.indexOf("chr")+5);
    //String pop1pop2 = name.substring(0,name.indexOf("_chr"));

    public static String listToString(List list, String separator) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < list.size(); i++) {
            if (i == list.size() - 1) {
                sb.append(list.get(i));
            } else {
                sb.append(list.get(i));
                sb.append(separator);
            }
        }
        return sb.toString();
    }

    /**
     * 将v1.0版本的基因中转换成 v1.1版本的gene名字，注意不是转录本
     * @param gene
     * @return
     */
    public static String getv11geneName(String gene){
        String out = null;
        if (gene.length() > 18){
            out = "NA";
            System.out.println("This gene is bad");
        }
        StringBuilder sb = new StringBuilder();
        Character goal = null;
        for (int i = 0; i < gene.length(); i++) {
            if (i==10){
                goal = '2';
                sb.append(goal);
            }
            else if (!(i==10)){
                goal = gene.charAt(i);
                sb.append(goal);
            }
        }
        out = sb.toString();
        return out;
    }

    /**
     * 将v1.1版本的基因中转换成 v1.0版本的gene名字，注意不是转录本
     * @param gene
     * @return
     */
    public static String getv10geneName(String gene){
        String out = null;
        if (gene.length() > 18){
            out = "NA";
            System.out.println("This gene is bad");
        }
        StringBuilder sb = new StringBuilder();
        Character goal = gene.charAt(10);
        for (int i = 0; i < gene.length(); i++) {
            if (i==10){
                goal = '1';
                sb.append(goal);
            }
            else if (!(i==10)){
                goal = gene.charAt(i);
                sb.append(goal);
            }
        }
        out = sb.toString();
        return out;
    }


}
