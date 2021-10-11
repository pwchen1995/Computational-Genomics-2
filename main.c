#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX 5

struct nodes{
    int id; // suffix id for leaf
    int parent_edge_label_start; // start and end
    int parent_edge_label_end;
    // child ptr
    struct nodes *child[MAX];
    struct nodes *parent;
    struct nodes *SL;
};

typedef struct nodes NODE;

NODE *root = NULL; // global variable
NODE *u = NULL, *v = NULL, *up = NULL, *vp = NULL;
NODE *temp = NULL,*temp2 = NULL,*temp3 = NULL;
int beta_remain = 0, beta = 0, internal_node = 0,child_c = 0,string_length = 0;
int BWTindex[100],a[11];
char alphabet2[] = {'$','A','C','G','T'};
char alphabet3[] = {'$','A','B','N'};
char alphabet1[] = {'$','I','M','P','S'};


NODE *creat(int,int,int);//給定節點資料內容，產生節點副函式
int split_info(char *[]);
char * split_seq(char *[]);
void insert_suffix(int,int,char[]);
void find_path(NODE *,int,int,int,char[]);
void NodeHop(NODE *,int,int,int,int,char[]);
void suffix_link(int,int,char[]);
void BWT_index(NODE *,char[]);
void BWT_loop(NODE *);
int c_compare(int ,char[],char[]);

int main(){

    FILE *input,*output;

    input = fopen("/Users/kevinchen/Desktop/WSU_Spring2020/Computational Genomics/project2/input3.fasta","r");
    output = fopen("/Users/kevinchen/Desktop/WSU_Spring2020/Computational Genomics/project2/BWT.txt","wb");
    int i = 0, j = 0, seq_temp_count = 0,info_length = 0;
    char *seq_whole[30000], *seq_info, *seq, bwt_index[11];
    // read file
    if(input != NULL) {
        while (fscanf(input, "%c", &seq_whole[i]) != EOF) {
            i++;
        }
    }
    //extract sequence
    info_length = split_info(seq_whole);
    seq = split_seq(seq_whole);
    printf("Length: %d Sequence: %s\n",i-info_length-1,seq);




    string_length = i-info_length;
    char seq_alt[sizeof(seq)], seq_temp[i-2];
    i = 0;
    while(i < string_length) {
        printf("Iteration %d\n", i);
        for (j = i; j < string_length; j++) {
            if(seq[j] == '\n'){
                break;
            }
            strncpy(&seq_temp[seq_temp_count], &seq[j], 1);
            printf("%s", seq_temp);
        }
        printf("\n");
        suffix_link(i, string_length - 1, seq); //傳入 id end位置, 與全部陣列做比對
        seq_temp_count = 0;
        i++;
    }
    printf("Internal Nodes: %d(including root)\n",internal_node+1);
    printf("BWT index:\n");
    BWT_index(root,seq);
    for (int k = 0; k <= 11; k++) {
        if (seq[a[k] - 1] == 0) {
            bwt_index[k] = seq[string_length - 1];
        } else {
            bwt_index[k] = seq[a[k] - 1];
        }
    }
    fprintf(output,"%s", bwt_index);
//    printf("%c",a);
//    fwrite(a, sizeof(int), sizeof(a), output);

    fclose(input);
    fclose(output);
//    system("PAUSE");
    return 0;
}


NODE *creat(int id, int start, int parentlabel){
    NODE *New; // 宣告一個節點指標
    New = (struct nodes *)malloc(sizeof(struct nodes));// creat new node
    New->id = id;
    New->parent_edge_label_start = start;
    New->parent_edge_label_end = parentlabel;
    New->child[0] = NULL;
    New->child[1] = NULL;
    New->child[2] = NULL;
    New->child[3] = NULL;
    New->child[4] = NULL;
    New->parent = NULL;
    New->SL = NULL;
    printf("Check S:%d E:%d\n",New->parent_edge_label_start,New->parent_edge_label_end);
    return (New); // 傳回New 指標
}
int split_info(char *seq_whole[]){
    int i = 0, count = 0, status = 0,seq_info_count = 0, seq_count = 0;
    char seq_info[50000];
    while(seq_whole[count] != '\0') {
        strncpy(&seq_info[seq_info_count], &seq_whole[count], 1);
        seq_info_count++;
        count++;
        if (seq_whole[count] == '\n') { count++;break; }
    }
    return  count;
}
char * split_seq(char *seq_whole[]){
    int count = 0, status = 0, seq_count = 0;
    char seq[50000];
    while(seq_whole[count] != '\0'){
        if(status == 0){
            if(seq_whole[count] == '\n'){
                status = 1;
            } // set status to 1 for sequence
        }
        if(seq_whole[count] != '\n' && status == 1){
            strncpy(&seq[seq_count],&seq_whole[count],1);
            seq_count++;
        }
        count++;
    }
    strncpy(&seq[seq_count],"$",1);
    return  seq;
}
void suffix_link(int suffid ,int parentlabel, char seq[]){
    //siffid = starting point , parentlabel = ending
    int cases = 0;
    int start = suffid, end = parentlabel; //assign original en and start to *start and *end

    if(root == NULL) { //creat root node
        root = creat(0,0,0);
        root -> SL = root;
        printf("Creat root! Add: %p\n",root);
        printf("----------------------------------------------------------------------------------\n");
        u = root;
        find_path(root,suffid,start,end,seq);
    }
    else{
        up = u->parent;
        if(u == root && u -> SL != NULL){cases = 2;} //case 1B
        else if (u != root && u -> SL != NULL){cases = 1;} // case 1A
        else if (up != root && u -> SL == NULL){cases = 3;} // case 2A
        else if (up == root && u -> SL == NULL){cases = 4;} // case 2B
        switch(cases){
            case 1: // case 1A
                v = u -> SL;
                find_path(v,suffid,start,end,seq);
                break;
            case 2: // case 1B
                v = u -> SL;
                find_path(v,suffid,start,end,seq);
                break;
            case 3: // case 2A
                vp = up -> SL;
                beta = u->parent_edge_label_end - u->parent_edge_label_start + 1; //beta
//                beta_start = u->parent_edge_label_start;
                NodeHop(vp,suffid,start,end,beta,seq);
//                find_path(v,suffid,start,end,seq);
                break;
            case 4: // case 2B
                vp = up -> SL;
                beta = u->parent_edge_label_end - u->parent_edge_label_start; //beta beta prime
//                beta_start = u->parent_edge_label_start;
                NodeHop(root,suffid,start,end,beta,seq);
//                find_path(v,suffid,start,end,seq);
                break;
            default:
                printf("Error input!\n");
                break;
        }
    }
}
void find_path(NODE * node, int suffid, int parentlabel_start, int parentlabel, char seq[]){
    NODE *ptr;
    int count_similar = 0, c = 0 , str_ptr = 0;
    int ori_start = 0,o = 0, s = 0,oo = 0, ss = 0 ,insert_pos = 0;//o = ascii of existing, s = current suffix
    int start = parentlabel_start, end = parentlabel,length_of_node = 0,length_of_suffix = 0,suffix_start = 0; //assign original en and start to *start and *end

    ptr = node;

    if (ptr->child[0] == NULL && ptr->child[1] == NULL&& ptr->child[2] == NULL&& ptr->child[3] == NULL&& ptr->child[4] == NULL) { //first insertion
        c = c_compare(start,alphabet2,seq);
        temp = creat(suffid, start, end);
        ptr->child[c] = temp;
        temp->parent = ptr;
        printf("Insert success! ID: %d, Start: %d ", temp->id, temp->parent_edge_label_start);
        printf("End: %d Add: %p Parent: %p Sufflink: %p\n", temp->parent_edge_label_end, temp, temp->parent,temp->SL);
        printf("---------------------------------------------------------------------------------------------------\n");
        u = temp -> parent;
    }
    else{ //
        c = c_compare(start,alphabet2,seq);
        if(ptr->child[c] == NULL){ // 要插入的child是空的
            temp = creat(suffid,start,end);
            ptr->child[c] = temp;
            temp ->parent = ptr;
            printf("Insert success! ID: %d, Start: %d ", temp->id, temp->parent_edge_label_start);
            printf("End: %d Add: %p Parent: %p Sufflink: %p\n", temp->parent_edge_label_end, temp, temp->parent,temp->SL);
            printf("---------------------------------------------------------------------------------------------------\n");
            u = temp->parent;
        }
        else { // child != NULL 1. right path 2. insert here
            //移動到下一個node找長度，如果長度大於suffix長度，findpath insert new node，如果於suffix長度大於node長度從這個node findpath
            ptr = ptr->child[c];
            length_of_node = ptr->parent_edge_label_end - ptr->parent_edge_label_start + 1;
            length_of_suffix = end - start + 1;
            if(length_of_node > length_of_suffix){
                ori_start = ptr->parent_edge_label_start; //pointer of existing child
                suffix_start = start; //pointer of suffix
                while(1){
                    o = seq[ori_start];
                    s = seq[suffix_start]; //walk down
                    if(o == s){
                        count_similar++;
                        ori_start++;
                        suffix_start++;
                    }
                    else{
                        temp2 = ptr; //hold original child
                        temp = creat(suffid,ptr->parent_edge_label_start,ptr->parent_edge_label_start+count_similar-1);
                        ptr = ptr->parent; //move back to parent and continue linking internal node
                        ptr->child[c] = temp;
                        temp->parent = ptr;

                        printf("Internal node Insert success! ID: %d, Start: %d ", temp->id, temp->parent_edge_label_start);
                        printf("End: %d Add: %p Parent: %p Sufflink: %p\n", temp->parent_edge_label_end, temp, temp->parent,temp->SL);
                        printf("---------------------------------------------------------------------------------------------------\n");

                        internal_node++; // this is the case of having internal node

                        c = c_compare(start+count_similar,alphabet2,seq);
                        temp3 = creat(suffid,start+count_similar,end);
                        temp->child[c] = temp3;
                        temp3->parent = temp;

                        printf("Link success! ID: %d, Start: %d ", temp3->id, temp3->parent_edge_label_start);
                        printf("End: %d Add: %p Parent: %p Sufflink: %p\n", temp3->parent_edge_label_end, temp3, temp3->parent,temp3->SL);
                        printf("---------------------------------------------------------------------------------------------------\n");

                        c = c_compare(ori_start,alphabet2,seq);
                        temp->child[c] = temp2;
                        temp2->parent = temp;

                        printf("Exist Link success! ID: %d, Start: %d ", temp2->id, temp2->parent_edge_label_start);
                        printf("End: %d Add: %p Parent: %p Sufflink: %p\n", temp2->parent_edge_label_end, temp2, temp2->parent,temp2->SL);
                        printf("---------------------------------------------------------------------------------------------------\n");

                        u = temp3->parent;
                        break;
                    }
                }
            }
            else{
//                c = c_compare(ptr->parent_edge_label_end+1,alphabet2,seq);
                start = start + length_of_node;
                find_path(ptr,suffid,start,end,seq);
            }
        }
    }
}
void NodeHop(NODE *node, int suffid, int parentlabel_start, int parentlabel, int beta,char seq[]) {
    NODE *ptr;
    int c = 0 ,insert_pos = 0;
    int start = parentlabel_start, end = parentlabel,length_of_node = 0,ori_start = 0,hop_start = 0; //assign original en and start to *start and *
    ptr = node;
    if(beta == 0){
        v = up;
        u->SL = v;
        find_path(v,suffid,parentlabel_start,parentlabel,seq);
    }
    else{
        c = c_compare(start,alphabet2,seq); //nodehop無論如何都有路徑
        ptr = ptr->child[c];
        length_of_node = ptr->parent_edge_label_end - ptr->parent_edge_label_start + 1;
//        start = start + length_of_node;
        hop_start = start + length_of_node;
        beta_remain = beta - length_of_node;
        if(length_of_node == beta) { //exist
            c = c_compare(start,alphabet2,seq);
            v = ptr;
            u->SL = v;
            find_path(v,suffid,start,end,seq);
        }
        else if(beta_remain > 0){ //need to creat
            c = c_compare(hop_start,alphabet2,seq);
            ptr = ptr->child[c];
            if(ptr->parent_edge_label_end-ptr->parent_edge_label_start+1 > beta_remain){
                ptr = ptr->parent;
                v = creat(suffid,ptr->parent_edge_label_end+1,ptr->parent_edge_label_end+beta_remain);
                u->SL = v;
                temp2 = ptr->child[c];
                ptr->child[c] = v;

                internal_node++; // this is the case of having internal node

                c = c_compare(temp2->parent_edge_label_start+beta,alphabet2,seq);
                v->child[c] = temp2;
                //parent 連回去
                v->parent = ptr;
                temp2->parent = v;
                find_path(v,suffid,start+beta,end,seq);
            }
            else{
                ptr = ptr->parent;
                NodeHop(ptr,suffid,start+ptr->parent_edge_label_end-ptr->parent_edge_label_start+1,end,beta_remain,seq);
            }
        }
        else{
            ori_start = ptr->parent_edge_label_start;
            v = creat(suffid,ori_start,ori_start+beta-1);
            u->SL = v;

            temp2 = ptr;
            ptr = ptr->parent;
            ptr->child[c] = v;
            v->parent = ptr;

            printf("Internal node Insert success! ID: %d, Start: %d ", v->id,v->parent_edge_label_start);
            printf("End: %d Add: %p Parent: %p Sufflink: %p\n", v->parent_edge_label_end, v, v->parent,v->SL);
            printf("---------------------------------------------------------------------------------------------------\n");

            internal_node++; // this is the case of having internal node

            c = c_compare(ori_start+beta,alphabet2,seq);
            v->child[c] = temp2;
            temp2->parent = v;

            printf("Exist Insert success! ID: %d, Start: %d ", temp2->id,temp2->parent_edge_label_start);
            printf("End: %d Add: %p Parent: %p Sufflink: %p\n", temp2->parent_edge_label_end, temp2, temp2->parent,temp2->SL);
            printf("---------------------------------------------------------------------------------------------------\n");

            u->SL = v;

            find_path(v,suffid,start+beta,end,seq);
        }
    }
}
void BWT_index(NODE *node,char seq[]) {
    NODE *ptr;
    int check = 0,root_c = 0,x = 0,check_c = 0;
    ptr = node;

    for (int i = 0; i < 5; i++) {
       if (ptr == NULL) {continue;}
       if (ptr->child[0] == NULL && ptr->child[1] == NULL && ptr->child[2] == NULL && ptr->child[3] == NULL && ptr->child[4] == NULL) {
                a[child_c] = ptr->id;
                child_c++;
                if(ptr->id == 0){
                    printf("%c ",seq[string_length-1]);
                }
                else{
                    printf("%c ",seq[ptr->id-1]);
                }
                i = 5;
       }
       else{
           if(ptr->child[i] == NULL){
               continue;
           }
           else{
               for (int j = 0; j < child_c; j++) {
                   if(ptr->child[i]->id== a[j]){check = 1;break;}
                   else{check = 0;}
               }
               if(check == 1){ continue;}
               else{
                   ptr = ptr->child[i];
                   BWT_index(ptr,seq);
                   ptr = ptr->parent;
               }
           }
       }
    }
}
int c_compare(int start,char alphabet2[],char seq[]){
    NODE *ptr;
    int o,s,c = 0;
    while(1){ //find right c
        o = alphabet2[c];
        s = seq[start];
        if(o != s){c++;}
        else{break;}
    }
    return c;
} //find the right child