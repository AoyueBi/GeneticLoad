package AoUtils;

import java.util.List;

public class AoString {
    public AoString(){

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
