classdef MMLP
    properties
        n_inputs % 输入层大小
        n_outputs % 输出层大小
        n_hidden_layers % 隐藏层数量
        n_hidden_units % 每个隐藏层中的神经元数量
        hidden_weights % 隐藏层权重
        output_weights % 输出层权重
        learning_rate % 学习率
        momentum % 动量参数
        prev_hidden_weight_delta % 上一次隐藏层权重更新的梯度
        prev_output_weight_delta % 上一次输出层权重更新的梯度
    end
    
    methods
        %Construction
        function obj = MMLP(n_inputs, n_outputs, n_hidden_layers, n_hidden_units, learning_rate, momentum)
            obj.n_inputs = n_inputs;
            obj.n_outputs = n_outputs;
            obj.n_hidden_layers = n_hidden_layers;
            obj.n_hidden_units = n_hidden_units;
            obj.learning_rate = learning_rate;
            obj.momentum = momentum;
            
            % 初始化隐藏层和输出层权重
            obj.hidden_weights = cell(1, n_hidden_layers);
            obj.output_weights = randn(n_hidden_units(end), n_outputs);
            
            for i = 1:n_hidden_layers
                if i == 1
                    input_size = n_inputs;
                else
                    input_size = n_hidden_units(i-1);
                end
                obj.hidden_weights{i} = randn(input_size, n_hidden_units(i));
            end
            
            % 初始化动量参数
            obj.prev_hidden_weight_delta = cell(1, n_hidden_layers);
            obj.prev_output_weight_delta = zeros(n_hidden_units(end), n_outputs);
            for i = 1:n_hidden_layers
                obj.prev_hidden_weight_delta{i} = zeros(size(obj.hidden_weights{i}));
            end
        end
        
        function [output, hidden_outputs] = forward(obj, input)
            % 前向传播计算输出和每个隐藏层的激活值
            hidden_outputs = cell(1, obj.n_hidden_layers);
            for i = 1:obj.n_hidden_layers
                if i == 1
                    hidden_input = input * obj.hidden_weights{i};
                else
                    hidden_input = hidden_outputs{i-1} * obj.hidden_weights{i};
                end
                hidden_outputs{i} = sigmoid(hidden_input);
            end
            
            output = sigmoid(hidden_outputs{end} * obj.output_weights);
        end
        
        function train(obj, inputs, targets)
            % 计算每个样本的输出和每个隐藏层的激活值
            [output, hidden_outputs] = obj.forward(inputs);
            
            % 计算输出层和每个隐藏层的误差信号
            output_error = (output - targets) .* sigmoid_derivative(output);
            hidden_errors = cell(1, obj.n_hidden_layers);
            hidden_errors{end} = output_error * obj.output_weights' .* sigmoid_derivative(hidden_outputs{end});
            for i = obj.n_hidden_layers-1:-1:1
                hidden_errors{i} = hidden_errors{i+1} * obj.hidden_weights{i+1}' .* sigmoid_derivative(hidden_outputs{i});
            end
            
            % 更新输出层权重
            obj.prev_output_weight_delta = obj.momentum * obj.prev_output_weight_delta - obj.learning_rate * hidden_outputs{end}' * output_error;
            obj.output_weights = obj.output_weights + obj.prev_output_weight_delta;
        
            % 更新隐藏层权重
            for i = obj.n_hidden_layers:-1:1
                if i == 1
                    hidden_input = inputs;
                else
                    hidden_input = hidden_outputs{i-1};
                end
                obj.prev_hidden_weight_delta{i} = obj.momentum * obj.prev_hidden_weight_delta{i} - obj.learning_rate * hidden_input' * hidden_errors{i};
                obj.hidden_weights{i} = obj.hidden_weights{i} + obj.prev_hidden_weight_delta{i};
            end
        end
    
        function output = predict(obj, inputs)
            % 预测给定输入的输出
            [output, ~] = obj.forward(inputs);
        end
    
        function set_learning_rate(obj, learning_rate)
            % 设置学习率
            obj.learning_rate = learning_rate;
        end
    
        function set_momentum(obj, momentum)
            % 设置动量参数
            obj.momentum = momentum;
        end
    end
end

function y = sigmoid(x)
    % sigmoid 激活函数
    y = 1./(1+exp(-x));
end

function y = sigmoid_derivative(x)
    % sigmoid 导数
    y = sigmoid(x) .* (1-sigmoid(x));
end
