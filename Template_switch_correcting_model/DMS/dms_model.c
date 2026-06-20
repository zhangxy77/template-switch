#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define MAX_ITER 10000
#define EPS 1e-10

// 辅助函数：计算对数似然
double calculate_likelihood(double *params, int K, double *S_matrix, double *Y_obs, 
                            int cycle, double D, double *N_i, double *r_i, 
                            double *p_ij, double *q_i, double *q_i_other, 
                            double *X_i, double *Y_ij) {
    
    // 1. 参数拆解
    double alpha = params[K];
    for (int i = 0; i < K; i++) {
        N_i[i] = params[i];
    }
    
    // 2. 相对丰度 r_i = N_i / sum(N_i)
    double sum_N = 0.0;
    for (int i = 0; i < K; i++) {
        sum_N += N_i[i];
    }
    if (sum_N <= 0) sum_N = 1e-8;
    for (int i = 0; i < K; i++) {
        r_i[i] = N_i[i] / sum_N;
    }
    
    // 3. 初始化数组
    memset(p_ij, 0, K * K * sizeof(double));
    memset(q_i, 0, K * sizeof(double));
    memset(q_i_other, 0, K * sizeof(double));
    memset(X_i, 0, K * sizeof(double));
    memset(Y_ij, 0, K * K * sizeof(double));
    
    // 4. TS概率 p_ij = alpha * S_ij * r_j
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            p_ij[i * K + j] = alpha * S_matrix[i * K + j] * r_i[j];
        }
        // 计算 q_i 和 q_i_other
        double sum_p = 0.0;
        for (int j = 0; j < K; j++) {
            sum_p += p_ij[i * K + j];
        }
        q_i[i] = sum_p;
        q_i_other[i] = q_i[i] - p_ij[i * K + i];
        
        // 防止数值问题
        if (q_i_other[i] > 0.99) q_i_other[i] = 0.99;
        if (q_i_other[i] < 1e-8) q_i_other[i] = 1e-8;
    }
    
    // 5. 非嵌合分子数 X_i = N_i * (2 - q_i_other)^cycle
    for (int i = 0; i < K; i++) {
        X_i[i] = N_i[i] * pow(2.0 - q_i_other[i], cycle);
    }
    
    // 6. 嵌合分子数 Y_ij (i != j)
    for (int i = 0; i < K; i++) {
        if (q_i_other[i] > 1e-8) {
            double term_base = pow(2.0, cycle) - pow(2.0 - q_i_other[i], cycle);
            for (int j = 0; j < K; j++) {
                Y_ij[i * K + j] = (p_ij[i * K + j] * N_i[i] / q_i_other[i]) * term_base;
            }
        }
    }
    
    // 7. 总分子数 M_total = sum(X_i) + sum_{i!=j} Y_ij
    double M_total = 0.0;
    for (int i = 0; i < K; i++) {
        M_total += X_i[i];
    }
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            if (i != j) {
                M_total += Y_ij[i * K + j];
            }
        }
    }
    
    // 避免除零
    if (M_total <= 0) M_total = 1e-8;
    
    // 8. 计算log-likelihood（只计算i!=j的嵌合分子）
    double log_lik = 0.0;
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            if (i != j) {
                double expected = D * Y_ij[i * K + j] / M_total;
                double observed = Y_obs[i * K + j];
                
                if (expected > 0) {
                    if (observed > 0) {
                        log_lik += observed * log(expected) - expected;
                    } else {
                        log_lik += -expected;
                    }
                }
            }
        }
    }
    
    // 检查是否为有限值
    if (!isfinite(log_lik)) log_lik = -1e10;
    
    return log_lik;
}

// 目标函数包装器
typedef struct {
    int K;
    double *S_matrix;
    double *Y_obs;
    int cycle;
    double D;
    double *work_buffer1;  // N_i
    double *work_buffer2;  // r_i
    double *work_buffer3;  // p_ij (K*K)
    double *work_buffer4;  // q_i
    double *work_buffer5;  // q_i_other
    double *work_buffer6;  // X_i
    double *work_buffer7;  // Y_ij (K*K)
} likelihood_data;

double objective_function(double *params, void *data) {
    likelihood_data *ld = (likelihood_data*)data;
    return calculate_likelihood(params, ld->K, ld->S_matrix, ld->Y_obs, 
                                ld->cycle, ld->D, ld->work_buffer1, ld->work_buffer2,
                                ld->work_buffer3, ld->work_buffer4, ld->work_buffer5,
                                ld->work_buffer6, ld->work_buffer7);
}

// 计算数值梯度
void compute_numerical_gradient(double *params, double *grad, int n, 
                                double (*func)(double*, void*), void *data, double eps) {
    double *params_plus = (double*)malloc(n * sizeof(double));
    double f0 = func(params, data);
    
    for (int i = 0; i < n; i++) {
        memcpy(params_plus, params, n * sizeof(double));
        params_plus[i] += eps;
        double f1 = func(params_plus, data);
        grad[i] = (f1 - f0) / eps;
    }
    free(params_plus);
}

// 简单的梯度下降优化（带边界约束）
int optimize_parameters(double *params, int n, double (*func)(double*, void*), void *data,
                        double *lower, double *upper, double *opt_value, int max_iter) {
    
    double *grad = (double*)malloc(n * sizeof(double));
    double *x_new = (double*)malloc(n * sizeof(double));
    double step_size = 0.001;
    double f_old = func(params, data);
    int iter;
    
    printf("  开始优化...\n");
    printf("  初始对数似然: %.6f\n", f_old);
    
    for (iter = 0; iter < max_iter; iter++) {
        // 计算梯度
        compute_numerical_gradient(params, grad, n, func, data, 1e-6);
        
        // 更新参数（梯度上升，因为我们要最大化似然）
        double max_change = 0.0;
        for (int i = 0; i < n; i++) {
            x_new[i] = params[i] + step_size * grad[i];  // 正梯度=上升
            
            // 边界投影
            if (x_new[i] < lower[i]) x_new[i] = lower[i];
            if (x_new[i] > upper[i]) x_new[i] = upper[i];
            
            double change = fabs(x_new[i] - params[i]);
            if (change > max_change) max_change = change;
        }
        
        // 更新参数
        for (int i = 0; i < n; i++) {
            params[i] = x_new[i];
        }
        
        double f_new = func(params, data);
        
        // 自适应步长
        if (f_new > f_old) {
            step_size *= 1.05;  // 增加步长
            f_old = f_new;
        } else {
            step_size *= 0.5;   // 减小步长
            // 回滚
            for (int i = 0; i < n; i++) {
                params[i] = params[i] - step_size * grad[i];
                if (params[i] < lower[i]) params[i] = lower[i];
                if (params[i] > upper[i]) params[i] = upper[i];
            }
        }
        
        if (iter % 100 == 0) {
            printf("  迭代 %d: 对数似然 = %.6f, 步长 = %.6f\n", iter, f_old, step_size);
        }
        
        // 检查收敛
        if (max_change < 1e-6 && iter > 10) {
            printf("  收敛于迭代 %d\n", iter);
            break;
        }
    }
    
    *opt_value = f_old;
    free(grad);
    free(x_new);
    
    return 0;
}

// 读取矩阵文件（第一行: K K，然后是矩阵数据）
double* read_matrix_file(const char *filename, int *K) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("错误：无法打开文件 %s\n", filename);
        return NULL;
    }
    
    int rows, cols;
    if (fscanf(file, "%d %d", &rows, &cols) != 2) {
        printf("错误：读取矩阵维度失败 %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    *K = rows;
    double *matrix = (double*)malloc(rows * cols * sizeof(double));
    if (!matrix) {
        printf("错误：内存分配失败\n");
        fclose(file);
        return NULL;
    }
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (fscanf(file, "%lf", &matrix[i * cols + j]) != 1) {
                printf("错误：读取矩阵数据失败 %s\n", filename);
                free(matrix);
                fclose(file);
                return NULL;
            }
        }
    }
    
    fclose(file);
    return matrix;
}

// 读取向量文件（第一行: K，然后是向量数据）
double* read_vector_file(const char *filename, int *K) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("错误：无法打开文件 %s\n", filename);
        return NULL;
    }
    
    if (fscanf(file, "%d", K) != 1) {
        printf("错误：读取向量维度失败 %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    double *vector = (double*)malloc(*K * sizeof(double));
    if (!vector) {
        printf("错误：内存分配失败\n");
        fclose(file);
        return NULL;
    }
    
    for (int i = 0; i < *K; i++) {
        if (fscanf(file, "%lf", &vector[i]) != 1) {
            printf("错误：读取向量数据失败 %s\n", filename);
            free(vector);
            fclose(file);
            return NULL;
        }
    }
    
    fclose(file);
    return vector;
}

int main() {
    // 从parameters.txt读取配置
    char y_obs_file1[256], y_obs_file2[256], s_matrix_file[256];
    char freq_obs_file1[256], freq_obs_file2[256];
    char output_file[256];
    int cycle, K;
    double D1, D2;
    
    FILE *param_file = fopen("parameters.txt", "r");
    if (!param_file) {
        printf("错误：无法打开parameters.txt\n");
        return 1;
    }
    
    // 读取文件路径
    if (fscanf(param_file, "%s", y_obs_file1) != 1) {
        printf("错误：读取 y_obs_file1 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%s", y_obs_file2) != 1) {
        printf("错误：读取 y_obs_file2 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%s", s_matrix_file) != 1) {
        printf("错误：读取 s_matrix_file 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%s", freq_obs_file1) != 1) {
        printf("错误：读取 freq_obs_file1 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%s", freq_obs_file2) != 1) {
        printf("错误：读取 freq_obs_file2 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%s", output_file) != 1) {
        printf("错误：读取 output_file 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%d", &cycle) != 1) {
        printf("错误：读取 cycle 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%lf", &D1) != 1) {
        printf("错误：读取 D1 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%lf", &D2) != 1) {
        printf("错误：读取 D2 失败\n");
        fclose(param_file);
        return 1;
    }
    if (fscanf(param_file, "%d", &K) != 1) {
        printf("错误：读取 K 失败\n");
        fclose(param_file);
        return 1;
    }
    
    fclose(param_file);
    
    printf("配置信息:\n");
    printf("  K = %d, cycle = %d, D1 = %.0f, D2 = %.0f\n", K, cycle, D1, D2);
    printf("  Y_obs_rep1: %s\n", y_obs_file1);
    printf("  Y_obs_rep2: %s\n", y_obs_file2);
    printf("  freq_obs_rep1: %s\n", freq_obs_file1);
    printf("  freq_obs_rep2: %s\n", freq_obs_file2);
    printf("\n");
    
    // 读取共享的S_matrix
    int S_K;
    double *S_matrix = read_matrix_file(s_matrix_file, &S_K);
    if (!S_matrix || S_K != K) {
        printf("错误：S_matrix维度不匹配\n");
        return 1;
    }
    
    // 用于存储两个rep的结果
    double **N_results = (double**)malloc(2 * sizeof(double*));
    double **r_results = (double**)malloc(2 * sizeof(double*));
    double alpha_results[2];
    double lik_results[2];
    int conv_results[2];
    
    for (int rep = 0; rep < 2; rep++) {
        N_results[rep] = (double*)malloc(K * sizeof(double));
        r_results[rep] = (double*)malloc(K * sizeof(double));
    }
    
    // 对每个rep进行处理
    for (int rep = 1; rep <= 2; rep++) {
        printf("\n========================================\n");
        printf("处理 REP %d ...\n", rep);
        printf("========================================\n");
        
        // 选择rep特定的文件
        char *y_obs_file = (rep == 1) ? y_obs_file1 : y_obs_file2;
        char *freq_obs_file = (rep == 1) ? freq_obs_file1 : freq_obs_file2;
        double D = (rep == 1) ? D1 : D2;
        
        // 读取Y_obs
        int Y_K;
        double *Y_obs = read_matrix_file(y_obs_file, &Y_K);
        if (!Y_obs || Y_K != K) {
            printf("错误：读取%s失败或维度不匹配\n", y_obs_file);
            continue;
        }
        
        // 读取freq_obs作为初始频率
        int freq_K;
        double *freq_initial = read_vector_file(freq_obs_file, &freq_K);
        if (!freq_initial || freq_K != K) {
            printf("错误：读取%s失败\n", freq_obs_file);
            free(Y_obs);
            continue;
        }
        
        printf("使用观测频率作为初始值:\n");
        for (int i = 0; i < K; i++) {
            printf("  freq_initial[%d] = %.6f\n", i+1, freq_initial[i]);
        }
        
        // 分配参数空间（N_i + alpha）
        double *params = (double*)malloc((K + 1) * sizeof(double));
        double *lower = (double*)malloc((K + 1) * sizeof(double));
        double *upper = (double*)malloc((K + 1) * sizeof(double));
        
        // 设置初始值 N_initial = N_total * freq_initial
        double N_total = 100000000.0;  // 总分子数，可根据需要调整
        for (int i = 0; i < K; i++) {
            params[i] = N_total * freq_initial[i];
            if (params[i] < 1e-6) params[i] = 1e-6;
            lower[i] = 1e-6;
            upper[i] = 1e8;
        }
        params[K] = 1e-4;  // alpha_initial
        lower[K] = 1e-8;
        upper[K] = 0.1;
        
        // 准备工作缓冲区
        likelihood_data ld;
        ld.K = K;
        ld.S_matrix = S_matrix;
        ld.Y_obs = Y_obs;
        ld.cycle = cycle;
        ld.D = D;
        
        ld.work_buffer1 = (double*)malloc(K * sizeof(double));
        ld.work_buffer2 = (double*)malloc(K * sizeof(double));
        ld.work_buffer3 = (double*)malloc(K * K * sizeof(double));
        ld.work_buffer4 = (double*)malloc(K * sizeof(double));
        ld.work_buffer5 = (double*)malloc(K * sizeof(double));
        ld.work_buffer6 = (double*)malloc(K * sizeof(double));
        ld.work_buffer7 = (double*)malloc(K * K * sizeof(double));
        
        if (!ld.work_buffer1 || !ld.work_buffer2 || !ld.work_buffer3 || 
            !ld.work_buffer4 || !ld.work_buffer5 || !ld.work_buffer6 || !ld.work_buffer7) {
            printf("错误：工作缓冲区分配失败\n");
            free(params); free(lower); free(upper); free(Y_obs); free(freq_initial);
            continue;
        }
        
        // 优化
        double final_lik;
        int convergence = optimize_parameters(params, K + 1, objective_function, &ld, 
                                              lower, upper, &final_lik, MAX_ITER);
        
        // 提取结果并存储
        double sum_N = 0.0;
        for (int i = 0; i < K; i++) {
            N_results[rep-1][i] = params[i];
            sum_N += params[i];
        }
        for (int i = 0; i < K; i++) {
            r_results[rep-1][i] = N_results[rep-1][i] / sum_N;
        }
        alpha_results[rep-1] = params[K];
        lik_results[rep-1] = final_lik;
        conv_results[rep-1] = convergence;
        
        printf("\nREP %d 完成:\n", rep);
        printf("  对数似然 = %.6f\n", final_lik);
        printf("  alpha = %.8f\n", alpha_results[rep-1]);
        
        // 清理rep特定的数据
        free(params);
        free(lower);
        free(upper);
        free(Y_obs);
        free(freq_initial);
        free(ld.work_buffer1);
        free(ld.work_buffer2);
        free(ld.work_buffer3);
        free(ld.work_buffer4);
        free(ld.work_buffer5);
        free(ld.work_buffer6);
        free(ld.work_buffer7);
    }
    
    // 以表格形式输出结果
    FILE *output = fopen(output_file, "w");
    if (!output) {
        printf("错误：无法创建输出文件 %s\n", output_file);
        free(S_matrix);
        for (int rep = 0; rep < 2; rep++) {
            free(N_results[rep]);
            free(r_results[rep]);
        }
        free(N_results);
        free(r_results);
        return 1;
    }
    
    // 写入参数信息
    fprintf(output, "Parameters:\n");
    fprintf(output, "K = %d\n", K);
    fprintf(output, "cycle = %d\n", cycle);
    fprintf(output, "D1 = %.0f\n", D1);
    fprintf(output, "D2 = %.0f\n", D2);
    fprintf(output, "alpha_rep1 = %.8f\n", alpha_results[0]);
    fprintf(output, "alpha_rep2 = %.8f\n", alpha_results[1]);
    fprintf(output, "log_likelihood_rep1 = %.6f\n", lik_results[0]);
    fprintf(output, "log_likelihood_rep2 = %.6f\n", lik_results[1]);
    fprintf(output, "convergence_rep1 = %d\n", conv_results[0]);
    fprintf(output, "convergence_rep2 = %d\n\n", conv_results[1]);
    
    // 写入表头
    fprintf(output, "Species\tN_rep1\tr_rep1\tN_rep2\tr_rep2\n");
    
    // 写入每个物种的数据
    for (int i = 0; i < K; i++) {
        fprintf(output, "%d\t%.6f\t%.8f\t%.6f\t%.8f\n", 
                i+1, 
                N_results[0][i], r_results[0][i],
                N_results[1][i], r_results[1][i]);
    }
    
    fclose(output);
    
    // 清理共享数据
    free(S_matrix);
    for (int rep = 0; rep < 2; rep++) {
        free(N_results[rep]);
        free(r_results[rep]);
    }
    free(N_results);
    free(r_results);
    
    printf("\n========================================\n");
    printf("所有处理完成！结果已保存到 %s\n", output_file);
    printf("========================================\n");
    
    return 0;
}