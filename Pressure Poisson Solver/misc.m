clear
close all
clc

I = rand(16,1).'

m_even = 0:2:15

for i = 1:length(m_even)
    I_new_even(i) = I(m_even(i) + 1);
end

I_new_even

m_odd = 1:2:15

for i = 1:length(m_odd)
    I_new_odd(i) = I(m_odd(i) + 1);
end
I
I_new_odd