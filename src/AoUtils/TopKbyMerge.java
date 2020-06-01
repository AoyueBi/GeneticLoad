package AoUtils;

/**
 * @author AoyueBi
 *
 */
public class TopKbyMerge{
    //数据描述：适用于无序单个数组，快排过程法利用快速排序的过程来求Top k。
    //实现步骤：根据快排规则，选择一个数为基准（代码是以最后一个数）分化数据，并记录基准数最后的落点下标，最后判断下标和k-1值大小（下标从0开始），不相等就继续朝k-1数量方向分化：
    //
    //下标小于k-1，对下标右侧（partion，end）继续二分；
    //
    //下标大于k-1，对下标左侧（first，partion）继续二分；
    //
    //直到k个数为top，但这k个数并没有顺序。

    int partion(int a[],int first,int end){ // 0, 9,6  int a[]={2,20,3,7,9,1,17,18,0,4};
        int i=first;
        int main=a[end];
        for(int j=first;j<end;j++){
            if(a[j]< main){
                int temp=a[j];
                a[j]=a[i];
                a[i]=temp;
                i++;
            }
        }
        a[end]=a[i];
        a[i]=main;
        return i;
    }
    void getTopKMinBySort(int a[],int first,int end,int k){
        if(first<end){
            int partionIndex=partion(a,first,end);
            if(partionIndex==k-1)return;
            else if(partionIndex>k-1)getTopKMinBySort(a,first,partionIndex-1,k);
            else getTopKMinBySort(a,partionIndex+1,end,k);
        }
    }
    public static void main(String []args){
        int a[]={2,20,3,7,9,1,17,18,0,4};
        int k=6;
        new TopKbyMerge().getTopKMinBySort(a,0,a.length-1,k);
        for(int i=0;i<k;i++){
            System.out.print(a[i]+" ");
        }
    }
}
