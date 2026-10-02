// Benchmark designs for placement. All are small, synthesizable RTL.

// ISCAS-85 c17: 6 NAND2 gates (used for exact-optimum comparisons)
module c17(input N1, N2, N3, N6, N7, output N22, N23);
  wire N10 = ~(N1 & N3);
  wire N11 = ~(N3 & N6);
  wire N16 = ~(N2 & N11);
  wire N19 = ~(N11 & N7);
  assign N22 = ~(N10 & N16);
  assign N23 = ~(N16 & N19);
endmodule

module adder8(input [7:0] a, b, input cin, output [7:0] s, output cout);
  assign {cout, s} = a + b + cin;
endmodule

module adder32(input [31:0] a, b, output [32:0] s);
  assign s = a + b;
endmodule

module mult6(input [5:0] a, b, output [11:0] p);
  assign p = a * b;
endmodule

module mult8(input [7:0] a, b, output [15:0] p);
  assign p = a * b;
endmodule

module alu8(input [7:0] a, b, input [2:0] op, output reg [7:0] y, output zero);
  always @* case (op)
    3'd0: y = a + b;   3'd1: y = a - b;   3'd2: y = a & b;   3'd3: y = a | b;
    3'd4: y = a ^ b;   3'd5: y = a << 1;  3'd6: y = a >> 1;  default: y = (a < b) ? 8'd1 : 8'd0;
  endcase
  assign zero = (y == 8'd0);
endmodule

// 16-bit CRC (CCITT) over an 8-bit input per clock: sequential, with flops
module crc16(input clk, input rst, input [7:0] d, output reg [15:0] crc);
  integer i;
  reg [15:0] c;
  always @* begin
    c = crc;
    for (i = 7; i >= 0; i = i - 1)
      c = {c[14:0], 1'b0} ^ ((c[15] ^ d[i]) ? 16'h1021 : 16'h0000);
  end
  always @(posedge clk) crc <= rst ? 16'hFFFF : c;
endmodule

module counter_cmp(input clk, rst, input [11:0] limit, output reg [11:0] q, output hit);
  assign hit = (q == limit);
  always @(posedge clk) q <= (rst | hit) ? 12'd0 : q + 12'd1;
endmodule

module mult12(input [11:0] a, b, output [23:0] p);
  assign p = a * b;
endmodule

module mult16(input [15:0] a, b, output [31:0] p);
  assign p = a * b;
endmodule

module alu32(input [31:0] a, b, input [2:0] op, output reg [31:0] y, output zero);
  always @* case (op)
    3'd0: y = a + b;   3'd1: y = a - b;   3'd2: y = a & b;   3'd3: y = a | b;
    3'd4: y = a ^ b;   3'd5: y = a << b[4:0];  3'd6: y = a >> b[4:0];  default: y = {31'd0, a < b};
  endcase
  assign zero = (y == 32'd0);
endmodule
