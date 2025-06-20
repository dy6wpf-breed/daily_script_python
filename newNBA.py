import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import time
from datetime import datetime, timedelta
import tkinter as tk

# 配置信息
API_KEY = "b4f35f5b-ad2c-4841-8fc9-f78549cbc8db"
HARDEN_ID = 279  # 哈登的 player_id

# 配置重试策略
retry_strategy = Retry(
    total=3,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504]
)
adapter = HTTPAdapter(max_retries=retry_strategy)
http = requests.Session()
http.mount("https://", adapter)
http.mount("http://", adapter)

# NBA球队英文名称到中文名称的映射
TEAM_NAME_MAPPING = {
    "LA Clippers": "快船",
    "Los Angeles Clippers": "快船",
    "Los Angeles Lakers": "湖人",
    "Golden State Warriors": "勇士",
    "Phoenix Suns": "太阳",
    "Sacramento Kings": "国王",
    "Dallas Mavericks": "独行侠",
    "Houston Rockets": "火箭",
    "Memphis Grizzlies": "灰熊",
    "New Orleans Pelicans": "鹈鹕",
    "San Antonio Spurs": "马刺",
    "Denver Nuggets": "掘金",
    "Minnesota Timberwolves": "森林狼",
    "Oklahoma City Thunder": "雷霆",
    "Portland Trail Blazers": "开拓者",
    "Utah Jazz": "爵士",
    "Boston Celtics": "凯尔特人",
    "Brooklyn Nets": "篮网",
    "New York Knicks": "尼克斯",
    "Philadelphia 76ers": "76人",
    "Toronto Raptors": "猛龙",
    "Chicago Bulls": "公牛",
    "Cleveland Cavaliers": "骑士",
    "Detroit Pistons": "活塞",
    "Indiana Pacers": "步行者",
    "Milwaukee Bucks": "雄鹿",
    "Atlanta Hawks": "老鹰",
    "Charlotte Hornets": "黄蜂",
    "Miami Heat": "热火",
    "Orlando Magic": "魔术",
    "Washington Wizards": "奇才"
}

def get_chinese_team_name(english_name):
    return TEAM_NAME_MAPPING.get(english_name, english_name)

def is_clippers_game(home_team, away_team):
    clippers_names = ["LA Clippers", "Los Angeles Clippers"]
    return any(name in (home_team, away_team) for name in clippers_names)

def get_live_games(target_date):
    url = "https://api.balldontlie.io/v1/games"
    headers = {"Authorization": API_KEY}
    
    params = {
        "start_date": target_date,
        "end_date": target_date,
        "per_page": 100
    }
    
    try:
        response = requests.get(url, headers=headers, params=params)
        if response.status_code == 200:
            return response.json()
        return None
    except Exception as e:
        print(f"获取比赛数据失败: {e}")
        return None

def get_harden_stats(game_id):
    url = "https://api.balldontlie.io/v1/stats"
    headers = {"Authorization": API_KEY} # 添加授权头
    params = {
        "game_ids[]": [game_id],
        "player_ids[]": [279],
        "per_page": 100
    }
    
    try:
        response = http.get(url, headers=headers, params=params) # 使用 http.get
        if response.status_code == 200:
            data = response.json()
            if data['data']:
                stats = data['data'][0]
                return {
                    'points': stats.get('pts', 0),
                    'assists': stats.get('ast', 0),
                    'rebounds': stats.get('reb', 0),
                    'steals': stats.get('stl', 0),
                    'blocks': stats.get('blk', 0),
                    'minutes': stats.get('min', '0')
                }
        else:
            pass # 保持静默，不打印API错误
        return None
    except Exception as e:
        print(f"获取哈登数据失败: {e}")
        return None

def format_game_message(game, is_clippers=False):
    home_team = get_chinese_team_name(game['home_team']['full_name'])
    away_team = get_chinese_team_name(game['visitor_team']['full_name'])
    home_score = game['home_team_score']
    visitor_score = game['visitor_team_score']
    status = game['status']
    period = game.get('period', 0)
    game_time = game.get('time', '').strip() # 获取 'time' 属性，并去除空白符

    # 临时调试打印：查看原始的 status 和 time 字段
    print(f"Game ID: {game['id']}, Raw Status: '{status}', Raw Time: '{game_time}'")

    message_parts = []
    message_parts.append(f"{home_team} vs {away_team}  {home_score} - {visitor_score}")

    display_status_info = ""

    if status.lower() == "final":
        display_status_info = "已结束"
    elif game_time and game_time != " ": # 如果 'time' 字段有实际内容（不是空白或"Final"）
        display_status_info = f"{status} (剩余时间: {game_time})"
    elif status and (":" in status and ("am" in status.lower() or "pm" in status.lower() or "et" in status.lower())): # 检查状态是否包含起始时间，例如 "7:00 pm ET"
        display_status_info = f"未开始 ({status})"
    else:
        display_status_info = status # 默认为原始状态

    message_parts.append(f"状态: {display_status_info}")
    if period and status.lower() != "final": # 只在比赛未结束时显示节数
        message_parts.append(f"节数: {period}")

    if is_clippers:
        harden_stats = get_harden_stats(game['id'])
        if harden_stats:
            message_parts.append("\n哈登数据：") # 添加空行和标题，使用 \n
            message_parts.append(f"得分: {harden_stats['points']} 助攻: {harden_stats['assists']} 篮板: {harden_stats['rebounds']}")
            message_parts.append(f"抢断: {harden_stats['steals']} 盖帽: {harden_stats['blocks']} 时间: {harden_stats['minutes']}")
    
    return "\n".join(message_parts) + "\n\n" # 额外空行用于分隔多场比赛，使用 \n

def update_display(text_widget):
    target_date = "2024-03-15" # 暂时修改为已知有比赛的日期
    
    games_data = get_live_games(target_date)
    
    display_message = "今天没有比赛数据。" # 默认消息
    
    if games_data and 'data' in games_data and games_data['data']:
        all_games_messages = []
        
        # 优先显示快船比赛
        clippers_games = [game for game in games_data['data'] if is_clippers_game(game['home_team']['full_name'], game['visitor_team']['full_name'])]
        if clippers_games:
            for game in clippers_games:
                all_games_messages.append(format_game_message(game, is_clippers=True))

        # 然后显示其他比赛
        other_games = [game for game in games_data['data'] if not is_clippers_game(game['home_team']['full_name'], game['visitor_team']['full_name'])]
        if other_games:
            for game in other_games:
                all_games_messages.append(format_game_message(game, is_clippers=False))

        if all_games_messages:
            display_message = "".join(all_games_messages).strip() # 拼接所有比赛信息并去除首尾空格
        else:
            display_message = "今天没有比赛。"

    # 检查滚动条是否在底部
    scroll_at_bottom = text_widget.yview()[1] == 1.0 # 如果第二个元素是1.0，表示已在底部

    text_widget.config(state=tk.NORMAL) # 允许编辑
    text_widget.insert(tk.END, display_message + "\n\n---\n\n") # 插入新消息并添加分隔线
    text_widget.update_idletasks() # 强制更新视图
    
    # 只有当滚动条之前就在底部时，才自动滚动到最新内容
    if scroll_at_bottom:
        text_widget.see(tk.END)

    text_widget.after(30000, lambda: update_display(text_widget)) # 每30秒更新一次

def main():
    root = tk.Tk()
    root.title("") # 将窗口标题设置为空
    root.geometry("250x350") # 设置窗口大小
    root.attributes("-topmost", True) # 置顶窗口
    # root.wm_attributes('-transparentcolor', 'SystemButtonFace') # 设置特定颜色为透明
    # root.config(bg="SystemButtonFace") # 设置窗口背景颜色为透明色
    root.resizable(True, True) # 允许窗口调节大小

    # 创建一个 Frame 来容纳 Text 控件和滚动条
    text_frame = tk.Frame(root, bg="#F0F0F0", borderwidth=0, highlightthickness=0)
    text_frame.pack(pady=10, padx=10, fill=tk.BOTH, expand=True)

    # 使用 tk.Text 来显示可滚动的比赛信息
    game_info_text = tk.Text(
        text_frame, # 将 Text 控件放置在 Frame 中
        font=("Microsoft YaHei UI", 12), # 调整字体大小
        fg="#333333", # 设置文本颜色为深灰色
        bg="#F0F0F0", # 设置背景颜色为浅灰色（保持不透明）
        selectbackground="#AAAAAA", # 选中背景色
        selectforeground="#FFFFFF", # 选中文字颜色
        borderwidth=0, # 移除边框
        highlightthickness=0, # 移除高亮边框
        wrap=tk.WORD, # 自动换行
        state=tk.NORMAL # 确保 Text 控件处于 NORMAL 状态以允许滚动
    )
    game_info_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    # 添加滚动条
    scrollbar = tk.Scrollbar(text_frame, command=game_info_text.yview) # 滚动条命令指向 Text 视图
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    game_info_text.config(yscrollcommand=scrollbar.set) # Text 控件的滚动命令指向滚动条

    # 绑定鼠标滚轮事件以确保滚动功能
    def on_mouse_wheel(event):
        if event.delta > 0:
            game_info_text.yview_scroll(-1, "units") # 向上滚动
        else:
            game_info_text.yview_scroll(1, "units") # 向下滚动

    game_info_text.bind("<MouseWheel>", on_mouse_wheel) # Windows
    game_info_text.bind("<Button-4>", lambda event: game_info_text.yview_scroll(-1, "units")) # Linux/macOS scroll up
    game_info_text.bind("<Button-5>", lambda event: game_info_text.yview_scroll(1, "units")) # Linux/macOS scroll down

    # 阻止键盘输入，保持只读特性
    game_info_text.bind("<Key>", lambda e: "break")

    # 显式设置焦点以确保鼠标滚轮事件被捕获
    game_info_text.focus_set()

    # 添加关闭按钮
    close_button = tk.Button(
        root,
        text="关闭",
        command=root.destroy,
        font=("Microsoft YaHei UI", 10),
        bg="#DDDDDD", # 按钮背景色
        fg="#333333", # 按钮文字颜色
        activebackground="#CCCCCC", # 按钮点击时背景色
        activeforeground="#333333", # 按钮点击时文字颜色
        relief=tk.FLAT # 扁平化按钮外观
    )
    close_button.pack(pady=(0, 10)) # 距离底部10像素

    # 修改 update_display 函数的调用方式，确保传入正确的 Text 控件
    def update_display_wrapper():
        update_display(game_info_text)

    update_display_wrapper() # 首次更新
    
    root.mainloop()

if __name__ == "__main__":
    main()