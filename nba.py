import requests
import time
from datetime import datetime, timedelta
import telebot

# 配置信息
API_KEY = "b4f35f5b-ad2c-4841-8fc9-f78549cbc8db"
TELEGRAM_BOT_TOKEN = "8185726760:AAFneLgJce5Adq3b2ClO35___eTSfslDaOk"
TELEGRAM_CHAT_ID = "6411482511"

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

# 添加这个函数定义
def is_clippers_game(home_team, away_team):
    """检查是否是快船的比赛"""
    clippers_names = ["LA Clippers", "Los Angeles Clippers"]
    return any(name in (home_team, away_team) for name in clippers_names)

# 初始化Telegram机器人
bot = telebot.TeleBot(TELEGRAM_BOT_TOKEN)

def get_chinese_team_name(english_name):
    return TEAM_NAME_MAPPING.get(english_name, english_name)

def get_live_games():
    url = "https://api.balldontlie.io/v1/games"
    headers = {
        "Authorization": f"{API_KEY}"
    }
    
    us_date = (datetime.now() - timedelta(days=1)).strftime('%Y-%m-%d')
    
    params = {
        "start_date": us_date,
        "end_date": us_date,
        "per_page": 100
    }
    
    try:
        response = requests.get(url, headers=headers, params=params)
        if response.status_code == 200:
            return response.json()
        else:
            return None
    except Exception as e:
        print(f"获取数据失败: {e}")
        return None

def get_harden_stats(game_id):
    url = f"https://api.balldontlie.io/v1/stats"
    headers = {
        "Authorization": f"{API_KEY}"
    }
    params = {
        "game_ids[]": [game_id],
        "player_ids[]": [192], # 哈登的player_id
        "per_page": 100
    }
    
    try:
        response = requests.get(url, headers=headers, params=params)
        if response.status_code == 200:
            data = response.json()
            if data['data']:
                stats = data['data'][0]
                return {
                    'points': stats['pts'],
                    'assists': stats['ast'],
                    'rebounds': stats['reb'],
                    'steals': stats['stl'],
                    'blocks': stats['blk'],
                    'minutes': stats['min']
                }
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
    
    if is_clippers:
        message = (
            "🏀 === 快船比赛 === 🏀\n"
            f"★ {home_team} vs {away_team} ★\n"
            f"★ 比分: {home_score} - {visitor_score} ★\n"
            f"★ 状态: {status} ★\n"
        )
        if period:
            message += f"★ 当前节数: {period} ★\n"
            
        # 获取哈登数据
        harden_stats = get_harden_stats(game['id'])
        if harden_stats:
            message += "\n🎯 哈登数据 🎯\n"
            message += f"得分: {harden_stats['points']}\n"
            message += f"助攻: {harden_stats['assists']}\n"
            message += f"篮板: {harden_stats['rebounds']}\n"
            message += f"抢断: {harden_stats['steals']}\n"
            message += f"盖帽: {harden_stats['blocks']}\n"
            message += f"上场时间: {harden_stats['minutes']}\n"
    else:
        message = (
            f"{home_team} vs {away_team}\n"
            f"比分: {home_score} - {visitor_score}\n"
            f"状态: {status}\n"
        )
        if period:
            message += f"当前节数: {period}\n"
    
    return message

def send_games_to_telegram(games_data):
    if not games_data or 'data' not in games_data:
        bot.send_message(TELEGRAM_CHAT_ID, "没有找到比赛数据")
        return
    
    games = games_data['data']
    if not games:
        bot.send_message(TELEGRAM_CHAT_ID, "今天没有比赛")
        return
    
    # 准备消息
    message = (
        f"🏀 NBA比赛信息更新 🏀\n"
        f"本地时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"美国时间: {(datetime.now() - timedelta(days=1)).strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"{'='*30}\n\n"
    )
    
    # 先处理快船比赛
    clippers_games = []
    other_games = []
    
    for game in games:
        home_team = game['home_team']['full_name']
        away_team = game['visitor_team']['full_name']
        if is_clippers_game(home_team, away_team):
            clippers_games.append(game)
        else:
            other_games.append(game)
    
    if clippers_games:
        message += "🔥 今天快船有比赛！ 🔥\n"
        for game in clippers_games:
            message += "\n" + format_game_message(game, is_clippers=True)
    
    if other_games:
        message += "\n其他比赛：\n"
        for game in other_games:
            message += "\n" + format_game_message(game)
    
    # 发送消息到Telegram
    try:
        bot.send_message(TELEGRAM_CHAT_ID, message)
    except Exception as e:
        print(f"发送到Telegram失败: {e}")

def main():
    update_interval = 300  # 5分钟更新一次
    
    # 发送启动消息
    bot.send_message(TELEGRAM_CHAT_ID, "NBA比赛信息推送服务已启动！")
    print("开始获取NBA比赛信息并发送到Telegram...")
    
    while True:
        try:
            games = get_live_games()
            if games:
                send_games_to_telegram(games)
            else:
                bot.send_message(TELEGRAM_CHAT_ID, "无法获取比赛信息，请检查API密钥是否正确")
            
            time.sleep(update_interval)
        except Exception as e:
            error_message = f"程序出错: {str(e)}\n5分钟后重试..."
            print(error_message)
            bot.send_message(TELEGRAM_CHAT_ID, error_message)
            time.sleep(300)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        bot.send_message(TELEGRAM_CHAT_ID, "程序已停止运行")
        print("\n程序已停止运行")
    except Exception as e:
        bot.send_message(TELEGRAM_CHAT_ID, f"程序异常退出: {e}")
        print(f"程序异常退出: {e}")