package AoUtils;

import java.util.Arrays;
import java.util.List;

public class AoString {
    public AoString(){

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
