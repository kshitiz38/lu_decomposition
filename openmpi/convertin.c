#include <stdio.h>

void doubleToRawString(double x, char *buffer) {
    const size_t bytesInDouble = 8;

    union {
        double value;
        unsigned char bytes[bytesInDouble];
    } u;

    u.value = x;

    unsigned char *input = u.bytes;
    char *output = buffer;

    int i;
    for(i = 0; i < bytesInDouble; ++i) {
        sprintf(output, "%02hhX", *input);

        ++input;
        output += 2;
    }
}

double rawStringToDouble(const char *input) {
    const size_t bytesInDouble = 8;

    union {
        double value;
        unsigned char bytes[bytesInDouble];
    } u;

    unsigned char *output = u.bytes;

    int i;
    for(i = 0; i < bytesInDouble; ++i) {
        sscanf(input, "%02hhX", output);

        input += 2;
        ++output;
    }

    return u.value;
}

int get_dim(FILE *fp, FILE *fp_w){
    char ch = getc(fp);
    int count = 1;
    while(ch != '\n') {
        if(ch == ' '){
            while(1){
                ch = getc(fp);
                if(ch=='\n'){
                    fprintf(fp_w, "%d\n", count);
                    return count;
                } else if (ch!=' '){
                    break;
                }
            }
            count++;
        }
        ch = getc(fp);
    }
    fprintf(fp_w, "%d\n", count);
    return count;
}

int main(){

    FILE *fp;
    FILE *fp_w;
    const size_t bytesInDouble = 8;

    fp_w = fopen("dims", "wb");
    fp = fopen("matrix", "rb");


    // double tmp;
    // int count = 0;
    // while(fscanf(fp, "%f", &tmp) == 1){
    //     count++;
    //     char ch = getc(fp);
    //     while(ch==' '){
    //         ch = getc(fp);
    //         if(ch=='\n') break;
    //     }
    //
    //
    // }
    int count = get_dim(fp, fp_w);
    printf("Dimension of Input Matrix is: %d\n", count);
    fclose(fp);
    fclose(fp_w);


    fp_w = fopen("pmatrix", "wb");
    fp = fopen("matrix", "rb");
    double val;
    if (fp != NULL)
    {
        while ( fscanf(fp, "%lf", &val) > 0 )
        {
            char hexstring[2*bytesInDouble + 1];
            doubleToRawString(val, hexstring);
            fprintf(fp_w, "%s", hexstring);
            double conv = rawStringToDouble(hexstring);
            // printf("%f %s %f\n", val, hexstring, conv);
        }
    }
    fclose(fp);
    fclose(fp_w);
}
