`timescale 1ns/1ns 
module MIPSpipeline24(clk, reset);
input clk, reset;
wire PC_write = 1, IF_Flush = 0, IFID_write = 1;
wire [31:0] PC,PCin, ID_PC4;
wire [31:0] PC4,PCbeq;////,ID_PC4,EX_PC4;
////wire [31:0] PCbeq,PC4beq;//,PCj,PC4bnej,PCjr; // PC signals in MUX
wire [31:0] Instruction, inst_reg,ID_Instruction; // Output of Instruction Memory
wire [5:0] ID_Opcode,ID_Function,EX_Function; //// Opcode, Function
// Extend
wire [15:0] ID_imm16; //// immediate in I type instruction
wire [31:0] EX_imm16; ////Im16_Ext,EX_Im16_Ext;
wire [31:0] sign_ext_out;
// regfile
wire [4:0] ID_rs,ID_rt,ID_rd,EX_rs,EX_rt,EX_rd,EX_WriteRegister,MEM_WriteRegister,WB_WriteRegister;
wire [31:0] WB_WriteDataToReg, ID_ReadData1, ID_ReadData2, EX_ReadData1, EX_ReadData2,MEM_WriteData;////,ReadData1Out,ReadData2Out;
wire [31:0] Mux_RS,Mux_RT;
// ALU
wire [31:0] Bus_A_ALU,Bus_B_ALU,Bus_B_forwarded;
wire [31:0] EX_ALUResult,MEM_ALUResult,WB_ALUResult;
wire ZeroFlag;////, OverflowFlag, CarryFlag, NegativeFlag,notZeroFlag;
wire [31:0] MEM_ReadDataOfMem,WB_ReadDataOfMem;////,WriteDataOfMe;
//Control signals 
wire RegDst,ALUSrc,MemtoReg,RegWrite,MemRead,MemWrite,Branch;//SignZero,JRControl;
wire PCSrc,HazMuxSel,compare_out;
wire [7:0] CON_SIG,OUT_CON;
wire ID_RegDst,ID_ALUSrc,ID_MemtoReg,ID_RegWrite,ID_MemRead,ID_MemWrite;////ID_Branch,ID_JRControl;
wire EX_RegDst,EX_ALUSrc,EX_MemtoReg,EX_RegWrite,EX_MemRead,EX_MemWrite;////EX_Branch,EX_JRControl;
wire MEM_MemtoReg,MEM_RegWrite,MEM_MemRead,MEM_MemWrite;
wire WB_MemtoReg,WB_RegWrite;
wire [4:0] Shift_Amt;
wire [1:0] ALUop,ID_ALUop,EX_ALUop;
wire [2:0] ALUControl;
//wire beqControl,notbeqControl;
wire [1:0] ForwardA,ForwardB;
wire ForwardC,ForwardD;
wire [31:0] shiftleft2_beq_out;//, addout; 
// PC Write Enable, IF/ID Write Enable
////wire PC_WriteEn,IFID_WriteEn;
assign CON_SIG [0] = RegDst; assign CON_SIG [1] = ALUSrc; assign CON_SIG [2] = MemtoReg;
assign CON_SIG [3] = RegWrite; assign CON_SIG [4] = MemRead; assign CON_SIG [5] = MemWrite;
assign CON_SIG [6] = ALUop [1]; assign CON_SIG [7] = ALUop [0];

assign ID_RegDst= OUT_CON [0]; assign ID_ALUSrc = OUT_CON [1]; assign ID_MemtoReg = OUT_CON [2];
assign ID_RegWrite = OUT_CON [3]; assign ID_MemRead = OUT_CON [4]; assign ID_MemWrite = OUT_CON [5];
assign ID_ALUop [1] = OUT_CON [6]; assign ID_ALUop [0] = OUT_CON [7];

assign ID_Opcode = ID_Instruction[31:26];
assign ID_Function = ID_Instruction[5:0];
assign ID_rs = ID_Instruction[25:21];
assign ID_rt = ID_Instruction[20:16];
assign ID_rd = ID_Instruction[15:11];
assign ID_imm16= ID_Instruction[15:0];

Clockgen c1 (clk,reset);
//IF Stage
pc    pc1(  clk , PCin , PC , reset , PC_write);
pcadder    pcadder11 (PC,PC4);
InstructionMem InstructionMem1(Instruction, PC);
instr_mux instruction_mux (Instruction,inst_reg,IF_Flush);
IFID pipeline_reg1 (clk,IFID_write,PC4,inst_reg,ID_Instruction,ID_PC4);
//ID Stage
HazardUnit HAZ (ID_rs,ID_rt,EX_rt,EX_MemRead,PC_write,IFID_write,HazMuxSel,IF_Flush,PCSrc);
andTOGetPcscr andgate(Branch,compare_out,PCSrc);
Control   hhh(RegDst,ALUSrc,MemtoReg,RegWrite,MemRead,MemWrite,Branch,ALUop,ID_Opcode);
Compare Comparator( Mux_RS, Mux_RT,  compare_out);
ForwardUnit_Branch FUB(MEM_WriteRegister,ID_rs,ID_rt,Branch,ForwardC,ForwardD);
comparator_input_rs_mux comp_mux1(ForwardC,ID_ReadData1,MEM_ALUResult,Mux_RS);
comparator_input_rt_mux comp_mux2(ForwardD,ID_ReadData2,MEM_ALUResult,Mux_RT);
RegFile  RegFile111  (ID_rs, ID_ReadData1 , ID_rt , ID_ReadData2 , clk , WB_RegWrite , WB_WriteRegister ,WB_WriteDataToReg );
control_mux  CON_MUX(CON_SIG,HazMuxSel,OUT_CON);
shift_left_2    shift_left_2111( shiftleft2_beq_out , sign_ext_out);
sign_extend    sign_extend1111( sign_ext_out,ID_imm16);
Addaddress PCadder(shiftleft2_beq_out,ID_PC4,PCbeq);
input_pc_mux IN_PC (PCSrc,PCbeq,PC4,PCin);
IDEX pipeline_reg2(clk,ID_MemtoReg,ID_RegWrite,ID_MemRead,ID_MemWrite,ID_RegDst,ID_ALUop,ID_ALUSrc,
ID_ReadData1,ID_ReadData2,sign_ext_out,ID_rs,ID_rt,ID_rd,ID_Function,EX_MemtoReg,EX_RegWrite,EX_MemRead,EX_MemWrite,EX_RegDst,EX_ALUop,
EX_ALUSrc,EX_ReadData1,EX_ReadData2,EX_imm16,EX_rs,EX_rt,EX_rd,EX_Function);
//EX Stage
ForwardUnit FU(MEM_WriteRegister,WB_WriteRegister,EX_rs,EX_rt,MEM_RegWrite,WB_MemtoReg,ForwardA,ForwardB);
alurs_mux ALUinput1mux(ForwardA,EX_ReadData1,WB_WriteDataToReg,MEM_ALUResult,Bus_A_ALU);
alurt_mux ALUinput2mux(ForwardB,EX_ReadData2,WB_WriteDataToReg,MEM_ALUResult,Bus_B_forwarded);
mux  muxalu (Bus_B_forwarded,EX_imm16,ALUSrc,Bus_B_ALU);
muxx mux225 (EX_rt,EX_rd,EX_RegDst,EX_WriteRegister);
ALUControl_Block    ALUControl_Block1111  ( ALUControl,Shift_Amt, EX_ALUop, EX_Function,EX_imm16);
alu   alu11( Bus_A_ALU, Bus_B_ALU, ALUControl,Shift_Amt,EX_ALUResult,ZeroFlag);
EXMEM pipeline_reg3(clk,EX_MemtoReg,EX_RegWrite,EX_MemRead,EX_MemWrite,EX_ALUResult,EX_WriteRegister,Bus_B_forwarded,
MEM_MemtoReg,MEM_RegWrite,MEM_MemRead,MEM_MemWrite,MEM_ALUResult,MEM_WriteRegister,MEM_WriteData);
//MEM Stage 
data_memory   data_memory1 (clk,MEM_ALUResult,MEM_WriteData,MEM_MemWrite,MEM_MemRead,MEM_ReadDataOfMem);
MEMWB pipeline_reg4(clk,MEM_MemtoReg,MEM_RegWrite,MEM_ReadDataOfMem,MEM_ALUResult,MEM_WriteRegister,
WB_MemtoReg,WB_RegWrite,WB_ReadDataOfMem,WB_ALUResult,WB_WriteRegister);
//WB Stage
mux  mux11 (WB_ALUResult,WB_ReadDataOfMem,WB_MemtoReg,WB_WriteDataToReg);
//// Add add11(addout,shiftleft2_beq_out ,PC4);/// to calculate the new address that beq will branch to
////and andbneControl(beqControl,Branch,ZeroFlag);
////mux  mux22(PC4,addout,beqControl,PCin);
endmodule 


module Clockgen (Clock,reset);
output reg Clock;
output reg reset;
initial
begin
Clock=0; reset=1;
end
always
begin
   reset<=0;
#1 Clock<=1;
#1 Clock<=0;
end
endmodule

module pc(  clock , pcin , pcout , reset , Pcwrite);
input clock , reset,Pcwrite;////Pcwrite for pipeline stall
input [31:0] pcin ;
output  reg [31:0]  pcout;
initial 
begin 
  pcout =0;
 end
always @(posedge clock)
begin
if (Pcwrite)
begin
    if(reset)
        pcout <= 0 ;
    else
        pcout <= pcin ;
end
end
endmodule


module RegFile(ra1, rd1 , ra2 , rd2 , clk , RegWrite , wa ,wd );
input[4:0] ra1;
output[31:0] rd1;
input[4:0] ra2;
output[31:0] rd2;
input clk, RegWrite;
//input werf ;
input[4:0] wa;
input[31:0] wd;
reg [31:0] registers[0:31];
initial 
begin 
registers[0]<=32'b00000000000000000000000000000000;// has to be 0
registers[1]<=32'b00000000000000000000000000000000;//0
registers[2]<=32'b00000000000000000000000000000011;//3
registers[3]<=32'b00000000000000000000000000000111;//7
registers[4]<=32'b00000000000000000000000000000111;//7
registers[5]<=32'b00000000000000000000000000000011;//3
registers[6]<=32'b00000000000000000000000000000010;//2
registers[7]<=32'b00000000000000000000000000000011;//3
registers[8]<=32'b00000000000000000000000000000001;//1
registers[9]<=32'b00000000000000000000000000000001;
registers[10]<=32'b00000000000000000000000000000000;
registers[11]<=32'b00000000000000000000000000000001;
registers[12]<=32'b00000000000000000000000000000001;
registers[13]<=32'b00000000000000000000000000000001;
registers[14]<=32'b00000000000000000000000000000001;
registers[15]<=32'b00000000000000000000000000000001;
registers[16]<=32'b00000000000000000000000000000001;
registers[17]<=32'b00000000000000000000000000000001;
registers[18]<=32'b00000000000000000000000000000001;
registers[19]<=32'b00000000000000000000000000000011;//3
registers[20]<=32'b00000000000000000000000000000010;//2
registers[21]<=32'b00000000000000000000000000000001;//1
registers[22]<=32'b00000000000000000000000000000001;
registers[23]<=32'b00000000000000000000000000000001;
registers[24]<=32'b00000000000000000000000000000001;
registers[25]<=32'b00000000000000000000000000000001;
registers[26]<=32'b00000000000000000000000000000001;
registers[27]<=32'b00000000000000000000000000000001;
registers[28]<=32'b00000000000000000000000000000001;
registers[29]<=32'b00000000000000000000000000000001;
registers[30]<=32'b00000000000000000000000000000001;
registers[31]<=32'b00000000000000000000000000000001;
end
assign rd1 = registers[ra1];
assign rd2 = registers[ra2];
always@ ( posedge clk )
    if (RegWrite)
        registers[wa] <= wd;
endmodule

 
module alu(       
      input          [31:0]     a,          //src1  
      input          [31:0]     b,          //src2  
      input          [2:0]     alu_control,     //function sel
      input	     [4:0]     ShiftAmount,  //for the instrution sll
      output     reg     [31:0]     result,          //result       
      output zero  
   );  
 always @(*)  
 begin   
      case(alu_control)  
      3'b000: result <= a + b; // add  
      3'b001: result <= a - b; // sub  
      3'b010: result <= a & b; // and  
      3'b011: result <= a | b; // or
      3'b100: result <= a<<ShiftAmount; //sll  
      default:result <= a + b; // add  
      endcase  
 end  
 assign zero = (result==32'd0) ? 1'b1: 1'b0;  
 endmodule  

 module data_memory  ( input clk,input[31:0]mem_access_addr,input[31:0] mem_write_data, input  mem_write_en, input mem_read, output [31:0] mem_read_data );  
      integer i;  
      reg [31:0] ram [0:1023];  
      wire [31: 0] ram_addr = mem_access_addr[31:0];  
      initial 
           begin  
           for(i=0;i<1024;i=i+1)  
                ram[i] <= i;  
            end 
      always @(posedge clk) begin  
           if (mem_write_en)  
                ram[ram_addr] <= mem_write_data;  
      end  
      assign mem_read_data = (mem_read==1'b1) ? ram[ram_addr]: 32'd0;   
 endmodule   

module ALUControl_Block( ALUControl,Shift_amt, ALUOp, Function,Sign_Extend);
output reg [2:0] ALUControl;
output reg [4:0] Shift_amt;
input [1:0] ALUOp;
input [5:0] Function;
input [31:0] Sign_Extend;
wire [7:0] ALUControlIn;
assign ALUControlIn = {ALUOp,Function};
always @(ALUControlIn)
casex (ALUControlIn)
 //8'b11xxxxxx: ALUControl=2'b01;  
 8'b00xxxxxx: ALUControl<=3'b000;//lw sw
 8'b01xxxxxx: ALUControl<=3'b001;// beq
 8'b10100000: ALUControl<=3'b000;//add
 8'b10100010: ALUControl<=3'b001;//sub
 8'b10100100: ALUControl<=3'b010;//and
 8'b10100101: ALUControl<=3'b011;//or
 8'b10000000: begin
	ALUControl<=3'b100;//sll
	Shift_amt <= Sign_Extend [10:6];
	end
 default: ALUControl=3'b000;
 endcase
endmodule

module sign_extend(sOut32,sIn16);
output [31:0] sOut32;
input [15:0] sIn16;
assign sOut32 = {{16{sIn16[15]}},sIn16};
endmodule
// Shift left 2 module 
module shift_left_2(Out32, In32);
output [31:0] Out32;
input [31:0] In32;

assign Out32 = {In32[29:0],2'b00};
endmodule

module zero_extend(zOut32,zIn16);
output [31:0] zOut32;
input [15:0] zIn16;
assign zOut32 = {{16{1'b0 }},zIn16};
endmodule
//---------------------------------------------------------------------------------------------------

module mux (in1,in2,sel,out1);

input [31:0] in1,in2;
output reg [31:0] out1;
input sel;

always@ (in1 or in2 or sel)

case (sel)
1'b0 : out1 <= in1 ;
1'b1 : out1 <= in2 ;
endcase 

endmodule

module muxx (in1,in2,sel,out1);

input [4:0] in1,in2;
output reg [4:0] out1;
input sel;

always@ (in1 or in2 or sel)

case (sel)
1'b0 : out1 <= in1 ; 
1'b1 : out1 <= in2 ;
endcase 

endmodule

module InstructionMem(instruction, address);

input [31:0] address;
output [31:0] instruction;
reg [31:0]instrmem[0:1023];
reg [31:0] temp;
assign instruction [31:0]=temp[31:0];

always @(address)
begin
 temp <=instrmem[address>>2];
end

initial
begin
$readmemb("Firsttestcase.txt", instrmem);
end
endmodule

module Control(
RegDst,
ALUSrc,
MemtoReg,
RegWrite,
MemRead,
MemWrite,
Branch,
ALUOp,
Opcode
);

output RegDst,ALUSrc,MemtoReg,RegWrite,MemRead,MemWrite,Branch;
output [1:0] ALUOp;
input [5:0] Opcode;
reg RegDst,ALUSrc,MemtoReg,RegWrite,MemRead,MemWrite,Branch;
reg [1:0] ALUOp;
always @(*)
casex (Opcode)
 6'b000000 : begin // R - type
     RegDst = 1'b1;
     ALUSrc = 1'b0;
     MemtoReg= 1'b0;
     RegWrite= 1'b1;
     MemRead = 1'b0;
     MemWrite= 1'b0;
     Branch = 1'b0;
     ALUOp = 2'b10;
     //Jump = 1'b0;
    // SignZero= 1'b0;
    end
 6'b100011 : begin // lw - load word
     RegDst = 1'b0;
     ALUSrc = 1'b1;
     MemtoReg= 1'b1;
     RegWrite= 1'b1;
     MemRead = 1'b1;
     MemWrite= 1'b0;
     Branch = 1'b0;
     ALUOp = 2'b00;
    // Jump = 1'b0;
     // SignZero= 1'b0; // sign extend
    end
 6'b101011 : begin // sw - store word
     RegDst = 1'bx;
     ALUSrc = 1'b1;
     MemtoReg= 1'bx;
     RegWrite= 1'b0;
     MemRead = 1'b0;
     MemWrite= 1'b1;
     Branch = 1'b0;
     ALUOp = 2'b00;
     //Jump = 1'b0;
     //SignZero= 1'b0;
    end


 6'b000100 : begin // beq- branch if equal
     RegDst = 1'b0;
     ALUSrc = 1'b0;
     MemtoReg= 1'b0;
     RegWrite= 1'b0;
     MemRead = 1'b0;
     MemWrite= 1'b0;
     Branch = 1'b1;
     ALUOp = 2'b01;
    // Jump = 1'b0;
    // SignZero= 1'b0; // sign extend
    end

/* 6'b000101 : begin // bne - branch if not equal
     RegDst = 1'b0;
     ALUSrc = 1'b0;
     MemtoReg= 1'b0;
     RegWrite= 1'b0;
     MemRead = 1'b0;
     MemWrite= 1'b0;
     Branch = 1'b1;
     ALUOp = 2'b01;
     Jump = 1'b0;
     SignZero= 1'b0; // sign extend
    end*/
 /*6'b001110 : begin // XORI - XOR immidiate
     RegDst = 1'b0;
     ALUSrc = 1'b1;
     MemtoReg= 1'b0;
     RegWrite= 1'b1;
     MemRead = 1'b0;
     MemWrite= 1'b0;
     Branch = 1'b0;
     ALUOp = 2'b11;
     Jump = 1'b0;
     SignZero= 1'b1; // zero extend
    end
 6'b000010 : begin // j - Jump
     RegDst = 1'b0;
     ALUSrc = 1'b0;
     MemtoReg= 1'b0;
     RegWrite= 1'b0;
     MemRead = 1'b0;
     MemWrite= 1'b0;
     Branch = 1'b0;
     ALUOp = 2'b00;
     Jump = 1'b1;
     SignZero= 1'b0;
    end*/
 default : begin 
     RegDst = 1'b0;
     ALUSrc = 1'b0;
     MemtoReg= 1'b0;
     RegWrite= 1'b0;
     MemRead = 1'b0;
     MemWrite= 1'b0;
     Branch = 1'b0;
     ALUOp = 2'b10;
     //Jump = 1'b0;
     //SignZero= 1'b0;
    end
 
endcase
endmodule
/*module zero_extend(zOut32,zIn16);
output [31:0] zOut32;
input [15:0] zIn16;
assign zOut32 = {{16{1'b0 }},zIn16};
endmodule*/


module Add(S,A,B);
output [31:0] S;
input [31:0] A,B;
assign S = A + B;
endmodule

//---------------------------------------------------------------------------------------------------
 module pcadder(pc,pcout);
input [31:0] pc;
output[31:0] pcout ;
assign  pcout=pc+4;
endmodule

////////////////////////////////////////////PIPELINE\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//IF/ID:
module IFID(clock,IFIDWrite,PC_Plus4,Inst,InstReg,PC_Plus4Reg); 
    input [31:0] PC_Plus4,Inst;//PC_Plus4 is the normal PC, no idea why this naming 
    input clock,IFIDWrite; 
    output reg [31:0] InstReg, PC_Plus4Reg;
     
    initial begin 
        InstReg = 0; 
        PC_Plus4Reg = 0; 
    end  
    always@(posedge clock) 
    begin ////yasmine made a mux for the flush
       // if(flush) 
       // begin 
       //    InstReg <= 0; 
       //    PC_Plus4Reg <=0; 
       // end 
         if(IFIDWrite)//else nothing is going to change 
        begin 
           InstReg <= Inst; 
           PC_Plus4Reg <= PC_Plus4; 
        end 
    end  
endmodule 

//ID/EX: 
module IDEX(clock,MemtoReg,RegWrite,MemRead,MemWrite,RegDst,ALUop,ALUSrc,DataA,DataB,imm_value,RegRs,RegRt,RegRd,Function,Memtoreg_reg,RegWrite_reg,
MemRead_reg,MemWrite_reg,RegDst_reg,ALUop_reg,ALUSrc_reg,DataAreg,DataBreg,imm_valuereg,RegRsreg,RegRtreg,RegRdreg,Functionreg); 
    input clock; 
    input MemtoReg,RegWrite,MemRead,MemWrite,RegDst,ALUSrc;////,Branch
    input [1:0] ALUop;
    input [4:0] RegRs,RegRt,RegRd; 
    input [31:0] DataA,DataB,imm_value;
    input [5:0] Function;
    output reg Memtoreg_reg,RegWrite_reg,MemRead_reg,MemWrite_reg,RegDst_reg,ALUSrc_reg;/////,Branch_reg
    output reg [1:0] ALUop_reg; 
    output reg [4:0] RegRsreg,RegRtreg,RegRdreg; 
    output reg [31:0] DataAreg,DataBreg,imm_valuereg; 
    output reg [5:0] Functionreg;
    initial begin
        Functionreg=0;           Memtoreg_reg = 0;        RegWrite_reg = 0;
        MemRead_reg = 0;         MemWrite_reg = 0;        RegDst_reg = 0; 
        ALUop_reg = 2'b00;       ALUSrc_reg = 0;   
        DataAreg = 0;            DataBreg = 0;            imm_valuereg = 0; 
        RegRsreg = 0;            RegRtreg = 0;            RegRdreg = 0; 
    end 
    always@(posedge clock) 
    begin 
	Functionreg<=Function;
	Memtoreg_reg <= MemtoReg;
        RegWrite_reg <= RegWrite;
        MemRead_reg <= MemRead;
        MemWrite_reg <= MemWrite; 
        RegDst_reg <= RegDst;
        ALUop_reg[1] <= ALUop[1];
	ALUop_reg[0] <= ALUop[0];
        ALUSrc_reg <= ALUSrc;  
        DataAreg <= DataA;
        DataBreg <= DataB; 
        imm_valuereg <= imm_value; 
        RegRsreg <= RegRs; 
        RegRtreg <= RegRt; 
        RegRdreg <= RegRd; 
    end 
endmodule 

//EX/MEM: 
module EXMEM(clock,MemtoReg,RegWrite,MemRead,MemWrite,ALUOut,RegRD,WriteDataIn,Memtoreg_reg,RegWrite_reg,
MemRead_reg,MemWrite_reg,ALUreg,RegRDreg,WriteDataOut); 
   input clock; 
   input MemtoReg,RegWrite,MemRead,MemWrite;////,Branch;
   input [4:0] RegRD; 
   input [31:0] ALUOut,WriteDataIn; 
   output reg Memtoreg_reg,RegWrite_reg,MemRead_reg,MemWrite_reg;////,Branch_reg;
   output reg [31:0] ALUreg,WriteDataOut; 
   output reg [4:0] RegRDreg; 
  
   initial begin 
      Memtoreg_reg = 0;         RegWrite_reg = 0;
      MemRead_reg = 0;         MemWrite_reg = 0;        
      ALUreg=0; 	       WriteDataOut=0; 		RegRDreg=0; 
   end 
    
    always@(posedge clock) 
    begin 
        Memtoreg_reg <= MemtoReg;
        RegWrite_reg <= RegWrite;
        MemRead_reg <= MemRead;
        MemWrite_reg <= MemWrite;
        ALUreg <= ALUOut; 
        RegRDreg <= RegRD; 
        WriteDataOut <= WriteDataIn; 
    end 
endmodule 

//MEM/WB: 
module MEMWB(clock,MemtoReg,RegWrite,Memout,ALUOut,RegRD,Memtoreg_reg,RegWrite_reg,Memreg,ALUreg,RegRDreg); 
   input clock; 
   input MemtoReg,RegWrite;
   input [4:0] RegRD; 
   input [31:0] Memout,ALUOut; 
   output reg Memtoreg_reg,RegWrite_reg; 
   output reg [31:0] Memreg,ALUreg; 
   output reg [4:0] RegRDreg; 
 
   initial begin 
      Memtoreg_reg = 0;         RegWrite_reg = 0;
      Memreg = 0;       ALUreg = 0;       RegRDreg = 0;    
   end 
    always@(posedge clock) 
    begin 
        Memtoreg_reg <= MemtoReg;
        RegWrite_reg <= RegWrite;
        Memreg <= Memout; 
        ALUreg <= ALUOut; 
        RegRDreg <= RegRD; 
    end 
endmodule

//forwording unit

 module ForwardUnit(MEMRegRd,WBRegRd,EXRegRs,EXRegRt, MEM_RegWrite, WB_MemtoReg, ForwardA, ForwardB); 
   input[4:0] MEMRegRd,WBRegRd,EXRegRs,EXRegRt;  
   input MEM_RegWrite, WB_MemtoReg; 
   output [1:0] ForwardA, ForwardB;
   reg [1:0] Forward_A, Forward_B;

   assign ForwardA = Forward_A; 		assign ForwardB = Forward_B; 
    
   //Forward A 
   always@(MEM_RegWrite or MEMRegRd or EXRegRs or WB_MemtoReg or WBRegRd) 
   begin 
      if((MEM_RegWrite)&&(MEMRegRd != 0)&&(MEMRegRd == EXRegRs)) //add add ///from EX/MEM to Alu
          Forward_A <= 2'b10;  ////ALU result
      else if((WB_MemtoReg)&&(WBRegRd != 0)&&(WBRegRd == EXRegRs)&&(MEMRegRd != EXRegRs) ) //lw stall add or add something not relevent add
          Forward_A <= 2'b01; ////output of wb stage
      else 
          Forward_A <= 2'b00; //// Rs single cycle
   end 
 
   //Forward B 
   always@(WB_MemtoReg or WBRegRd or EXRegRt or MEMRegRd or MEM_RegWrite) 
   begin 
      if((WB_MemtoReg)&&(WBRegRd != 0)&&(WBRegRd == EXRegRt)&&(MEMRegRd != EXRegRt) ) //lw stall add or add something not relevent add
          Forward_B <= 2'b01;/////////output of wb stage 
      else if((MEM_RegWrite)&&(MEMRegRd != 0)&&(MEMRegRd == EXRegRt)) //add add
          Forward_B <= 2'b10;/////////ALU result 
      else  
          Forward_B <= 2'b00; //// Rt single cycle 
   end 
 
 endmodule
//Hazard Detection Unit Module: 
module HazardUnit(IDRegRs,IDRegRt,EXRegRt,EXMemRead,PCWrite,IFIDWrite,HazMuxCon,if_flush,pcsrcs/*,Mcont*/); 
    
    input [4:0] IDRegRs,IDRegRt,EXRegRt; 
    input EXMemRead,pcsrcs;
    //input [2:0] Mcont ;//branch mem read memwrit ..control sig going to memory stage if branch  001:mem read 0mem write 0 branch 1
    output PCWrite, IFIDWrite, HazMuxCon,if_flush; 
    reg PC_Write, IFID_Write, Haz_MuxCon,if__flush;
    assign  PCWrite = PC_Write; 		assign IFIDWrite = IFID_Write; 
    assign HazMuxCon = Haz_MuxCon;		assign if_flush = if__flush;
    always@(IDRegRs,IDRegRt,EXRegRt,EXMemRead,pcsrcs/*,Mcont*/) 
     begin
    if(EXMemRead&&((EXRegRt == IDRegRs)||(EXRegRt == IDRegRt)))////lw add 
       begin//stall 
           PC_Write <= 0; 
           IFID_Write <= 0; 
           Haz_MuxCon <= 1; //kda hatal3 al signal bta3t alcontrol unit =0
       end
    else if (pcsrcs)//((pcsrcs)&&(Mcont))
    begin
     if__flush <=1;
    end 
    else 
       begin//no stall 
            PC_Write <= 1; 
            IFID_Write <= 1; 
            Haz_MuxCon <= 0; 
            if__flush <=0;   
       end 
 end
endmodule

//datatowrite ..output of mux (selector memtoreg) in wb stage//memaluout ..output of ex/mem reg (ALUresult)
//aluscra is the output of alurs_mux and is the input 1 to alu
module alurs_mux(ForwardA,Rs,datatowrite,MEMALUOut,ALUSrcA); 
input [1:0] ForwardA; 
input [31:0] Rs,datatowrite,MEMALUOut; 
output reg [31:0] ALUSrcA; 
 
always@(ForwardA,Rs,datatowrite,MEMALUOut) 
begin 
  case(ForwardA) 
      2'b00: 
        ALUSrcA <= Rs; 
      2'b01: 
        ALUSrcA <= datatowrite; 
      2'b10: 
        ALUSrcA <= MEMALUOut; 
  endcase 
end //OUTPUT IS INPUT TO RS
 endmodule

module alurt_mux(ForwardB,Rt,datatowrite,MEMALUOut,ALUSrcB); 
input [1:0] ForwardB; 
input [31:0] Rt,datatowrite,MEMALUOut; 
output reg [31:0] ALUSrcB;  
 
always@(ForwardB,Rt,datatowrite,MEMALUOut) 
begin 
  case(ForwardB) 
      2'b00: 
        ALUSrcB <= Rt; 
      2'b01: 
        ALUSrcB <= datatowrite; 
      2'b10: 
        ALUSrcB <= MEMALUOut; 
  endcase 
end 
//OUTPUT OF THIS MUX IS INPUT TO MUX (IMMEDIATE & RT ) 
endmodule 


module ForwardUnit_Branch(MEMRegRd,IDRegRs,IDRegRt,branch,ForwardC, ForwardD); 
   input[4:0] MEMRegRd,IDRegRs,IDRegRt;  
   input branch; 
   output reg ForwardC, ForwardD;  
    
   //Forward C 
   always@(branch or MEMRegRd or IDRegRs  ) 
   begin 
      if((branch)&&(MEMRegRd != 0)&&(MEMRegRd == IDRegRs)) 
         ForwardC <= 1; 
      else ForwardC <= 0;
    end
   //Forward D 
   always@(branch or MEMRegRd or IDRegRt ) 
   begin 
      if((branch)&&(MEMRegRd != 0)&&(MEMRegRd == IDRegRt)) 
         ForwardD <= 1; 
      else ForwardD <= 0; 
   end
 endmodule


//compare rs if equal rt
module Compare(
    input  [31:0] RS,
    input  [31:0] RT,
    output EQ);
  assign EQ  = ( RS== RT);
endmodule
//calculate address of the label of the branch 
module Addaddress(
    input  [31:0] shiftleft2result,
    input  [31:0] pcplusfour,
    output [31:0] labelAddr
    );

    assign labelAddr = (shiftleft2result + pcplusfour);

endmodule
module comparator_input_rs_mux(Forwardc,Rs,ex_mem_alu_out,rsToComarator); 
input Forwardc; 
input [31:0] Rs,ex_mem_alu_out; 
output reg [31:0] rsToComarator; 
 
always@(Forwardc,Rs,ex_mem_alu_out) 
begin 
  case(Forwardc) 
      0: 
       rsToComarator<= Rs; 
      1: 
       rsToComarator <= ex_mem_alu_out; 
  endcase 
end //OUTPUT IS INPUT comparator rs
 endmodule

module comparator_input_rt_mux(Forwardd,Rt,ex_mem_alu_out1,rsToComarator1); 
input Forwardd; 
input [31:0] Rt,ex_mem_alu_out1; 
output reg [31:0] rsToComarator1;  
 
always@(Forwardd,Rt,ex_mem_alu_out1) 
begin 
  case(Forwardd) 
      0: 
       rsToComarator1<= Rt; 
      1: 
       rsToComarator1 <= ex_mem_alu_out1; 
  endcase 
end //OUTPUT IS INPUT comparator  rt
 endmodule

module andTOGetPcscr (branchCon,compOut,pcscr);//compOut comparator output .. branch con branch control signal
input branchCon,compOut;
output pcscr;
assign pcscr=(branchCon&compOut);
endmodule
////yasmin kalakna
module instr_mux (instruction,instr_reg, IF_flush);
input [31:0]instruction;
input IF_flush;
output reg [31:0] instr_reg;

always @(instruction, IF_flush)
begin
case(IF_flush) 
      1: 
       instr_reg <= 32'b0; 
      0: 
      instr_reg <= instruction; 
  endcase  
end 
endmodule

module muxdestination (reg_dstn,id_ex_rt,id_ex_rd, regDst_control);//choose between rt &rd
input[31:0]id_ex_rt,id_ex_rd;
input regDst_control;
output reg [31:0] reg_dstn;

always @ (id_ex_rt,id_ex_rd, regDst_control)
begin
case(regDst_control) 
    0: 
       reg_dstn<=id_ex_rt; 
     1: 
     reg_dstn <= id_ex_rd; 
 endcase  
end 
endmodule

module input_pc_mux (pcSrc,branch_pc,pcPLUSfour,ToPc);
input[31:0]branch_pc,pcPLUSfour;
input pcSrc;
output [31:0] ToPc;
reg [31:0] ToPc;
always @ (pcSrc,branch_pc,pcPLUSfour)
begin
case(pcSrc) 
    0: 
       ToPc<=pcPLUSfour; 
     1: 
     ToPc <= branch_pc; 
 endcase  
end 
endmodule
module control_mux (CON_SIG,HAZ_SEL,OUT_CON);
input [7:0] CON_SIG;
input HAZ_SEL;
output reg [7:0] OUT_CON;
always @ (CON_SIG,HAZ_SEL)
case(HAZ_SEL)
0:OUT_CON<=CON_SIG;
1:OUT_CON<=0;
endcase
endmodule