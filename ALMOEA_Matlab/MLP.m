classdef MLP < handle
    properties(SetAccess = private)
        input_size   % size of input
        hidden_size  % size of hidden layer
        output_size  % size of output
        hidden_weights  % weights between input and hidden
        output_weights  % weights between hidden and output
    end
    
    methods
        % Constructor
        function obj = MLP(input_size, hidden_size, output_size)
            obj.input_size = input_size;
            obj.hidden_size = hidden_size;
            obj.output_size = output_size;
            % Initialization of weights
            obj.hidden_weights = randn(hidden_size, input_size);
            obj.output_weights = randn(output_size, hidden_size);
        end
        
        % 前向传播
        function [hidden_output, output] = forward(obj, input)
            % 计算隐藏层输出
            hidden_output = sigmoid(obj.hidden_weights * input);
            % 计算输出层输出
            output = sigmoid(obj.output_weights * hidden_output);
        end
        
        % 反向传播
        function [grad_hidden_weights, grad_output_weights] = backward(obj, input, target, hidden_output, output)
            % 计算输出误差
            output_error = (target - output) .* output .* (1 - output);
            % 计算隐藏层误差
            hidden_error = obj.output_weights' * output_error .* hidden_output .* (1 - hidden_output);
            % 计算输出层权重梯度
            grad_output_weights = output_error * hidden_output';
            % 计算隐藏层权重梯度
            grad_hidden_weights = hidden_error * input';
        end
        
        % 更新权重
        function obj = update(obj, grad_hidden_weights, grad_output_weights, learning_rate)
            % 更新隐藏层权重
            obj.hidden_weights = obj.hidden_weights + learning_rate * grad_hidden_weights;
            % 更新输出层权重
            obj.output_weights = obj.output_weights + learning_rate * grad_output_weights;
        end
    end
end

function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end