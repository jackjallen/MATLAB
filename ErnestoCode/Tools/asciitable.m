function ascii
% ASCII shows a chart of ascii values
%
%  32 - 47  :   ! " # $ % & ' ( ) * + , - . / 
%  48 - 57  : 0 1 2 3 4 5 6 7 8 9 
%  56 - 64  : : ; < = > ? @ 
%  65 - 90  : A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 
%  91 - 97  : [ \ ] ^ _ ` a 
%  98 - 122 : b c d e f g h i j k l m n o p q r s t u v w x y z 
% 123 - 126 : { | } ~ 
% 127 - 159 :   ?            ?  ? ?             ?   
% 160 - 189 :   ¡ ¢ £ ¤ ¥ ¦ § ¨ © ª « ¬ ­ ® ¯ ° ± ² ³ ´ µ ¶ · ¸ ¹ º » ¼ ½ 
% 190 - 219 : ¾ ¿ À Á Â Ã Ä Å Æ Ç È É Ê Ë Ì Í Î Ï Ð Ñ Ò Ó Ô Õ Ö × Ø Ù Ú Û 
% 220 - 249 : Ü Ý Þ ß à á â ã ä å æ ç è é ê ë ì í î ï ð ñ ò ó ô õ ö ÷ ø ù 
% 250 - 256 : ú û ü ý þ ÿ 
%

fprintf(' 32 - 47  : '); fprintf('%c ', char( 32:47 )); fprintf('\n');
fprintf(' 48 - 57  : '); fprintf('%c ', char( 48:57 )); fprintf('\n');
fprintf(' 56 - 64  : '); fprintf('%c ', char( 58:64 )); fprintf('\n');
fprintf(' 65 - 90  : '); fprintf('%c ', char( 65:90 )); fprintf('\n');
fprintf(' 91 - 97  : '); fprintf('%c ', char( 91:97 )); fprintf('\n');
fprintf(' 98 - 122 : '); fprintf('%c ', char( 98:122)); fprintf('\n');
fprintf('123 - 126 : '); fprintf('%c ', char(123:126)); fprintf('\n');
fprintf('127 - 159 : '); fprintf('%c ', char(127:159)); fprintf('\n');
fprintf('160 - 189 : '); fprintf('%c ', char(160:189)); fprintf('\n');
fprintf('190 - 219 : '); fprintf('%c ', char(190:219)); fprintf('\n');
fprintf('220 - 249 : '); fprintf('%c ', char(220:249)); fprintf('\n');
fprintf('250 - 255 : '); fprintf('%c ', char(250:255)); fprintf('\n');
